# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module with helpers for broadcasting results to samp clients
"""
import contextlib
import os
import tempfile

from astropy.samp import SAMPIntegratedClient


__all__ = [
    'find_client_id', 'send_table_to', 'send_product_to', 'send_spectrum_to',
    'send_image_to', 'accessible_table', 'connection']


def find_client_id(conn, name):
    """returns the SAMP id of the client with samp.name samp_name.

    This will raise a KeyError if the client is not on the hub.
    """
    for client_id in conn.get_registered_clients():
        if conn.get_metadata(client_id).get("samp.name") == name:
            return client_id
    raise KeyError(name)


def send_table_to(conn, table, client_name=None, name="data"):
    """
    sends astropy_table via SAMP.
    """
    with accessible_table(table) as url:
        message = {
            "samp.mtype": "table.load.votable",
            "samp.params": {
                "url": url,
                "name": name,
            },
        }
        if client_name is None:
            for client_id in conn.get_registered_clients():
                conn.call_and_wait(client_id, message, "10")
        else:
            client_id = find_client_id(conn, client_name)
            conn.call_and_wait(client_id, message, "10")


def send_product_to(conn, url, mtype, client_name=None, name="data"):
    """
    sends SAMP messages to load data.

    This is a helper for send_spectrum_to and send_image_to, which work exactly
    analogous to each other, except that the mtypes are different.

    If dest_client_id, this is a broadcast (and we don't wait for any
    responses). If dest_client_id is given, we wait for acknowledgement by the
    receiver.
    """
    message = {
        "samp.mtype": mtype,
        "samp.params": {
            "url": url,
            "name": name,
        },
    }
    if client_name is None:
        conn.notify_all(message)
    else:
        client_id = find_client_id(conn, client_name)
        conn.notify(client_id, message)


def send_spectrum_to(conn, url, client_name=None, name="data"):
    """
    asks a spectrum client to open a remote spectrum via SAMP.
    """
    send_product_to(
        conn, url, "spectrum.load.ssa-generic",
        client_name=client_name, name=name)


def send_image_to(conn, url, client_name=None, name="data"):
    """
    asks an image client to open a remote image via SAMP.
    """
    send_product_to(
        conn, url, "image.load.fits",
        client_name=client_name, name=name)


@contextlib.contextmanager
def accessible_table(table):
    """
    a context manager making astropy_table available under a (file)
    URL for the controlled section.
    """
    handle, f_name = tempfile.mkstemp(suffix=".xml")
    with open(handle, "w") as f:
        table.write(output=f, format="votable")
    try:
        yield "file://" + f_name
    finally:
        os.unlink(f_name)


@contextlib.contextmanager
def connection(
    client_name="pyvo client", description="A generic PyVO client", **kwargs
):
    """
    a context manager to give the controlled block a SAMP connection.
    The program will disconnect as the controlled block is exited.
    """
    client = SAMPIntegratedClient(
        name=client_name, description=description, **kwargs)
    client.connect()
    try:
        yield client
    finally:
        client.disconnect()
