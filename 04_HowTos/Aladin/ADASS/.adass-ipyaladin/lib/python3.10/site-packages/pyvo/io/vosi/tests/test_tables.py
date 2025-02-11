#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.io.vosi
"""
import contextlib
import io
import pytest

import pyvo.io.vosi as vosi
import pyvo.io.vosi.vodataservice as vs
from pyvo.io.vosi.exceptions import (
    W02, W03, W04, W05, W06, W07, W08, W09, W10, W11, W12, W13, W14, W37)
from pyvo.io.vosi.exceptions import E01, E02, E03, E06

from astropy.utils.data import get_pkg_data_filename


class TestTables:
    def test_all(self):
        tablesfile = vosi.parse_tables(
            get_pkg_data_filename("data/tables.xml"))
        table = next(tablesfile.iter_tables())

        assert table.name == "test.all"
        assert table.title == "Test table"
        assert table.description == "All test data in one table"
        assert table.utype == "utype"
        assert table.nrows == 30

        col = table.columns[0]
        fkc = table.foreignkeys[0]

        assert col.name == "id"
        assert col.description == "Primary key"
        assert col.unit == "unit"
        assert col.ucd == "meta.id;meta.main"
        assert col.utype == "utype"

        assert isinstance(col.datatype, vs.TAPType)
        assert str(col.datatype) == "<DataType arraysize=*>VARCHAR</DataType>"
        assert col.datatype.arraysize == "*"
        assert col.datatype.delim == ";"
        assert col.datatype.size == "42"
        assert col.datatype.content == "VARCHAR"

        assert "indexed" in col.flags
        assert "primary" in col.flags

        assert fkc.targettable == "test.foreigntable"
        assert fkc.fkcolumns[0].fromcolumn == "testkey"
        assert fkc.fkcolumns[0].targetcolumn == "testkey"
        assert fkc.description == "Test foreigner"
        assert fkc.utype == "utype"

    def _test_datatypes_votable(self, cols):
        assert cols[0].datatype.content == 'boolean'
        assert cols[1].datatype.content == 'bit'
        assert cols[2].datatype.content == 'unsignedByte'
        assert cols[3].datatype.content == 'short'
        assert cols[4].datatype.content == 'int'
        assert cols[5].datatype.content == 'long'
        assert cols[6].datatype.content == 'char'
        assert cols[7].datatype.content == 'unicodeChar'
        assert cols[8].datatype.content == 'float'
        assert cols[9].datatype.content == 'double'
        assert cols[10].datatype.content == 'floatComplex'
        assert cols[11].datatype.content == 'doubleComplex'

    def test_datatypes_votable(self):
        tablesfile = vosi.parse_tables(
            get_pkg_data_filename("data/tables/datatypes_votable.xml"))

        votable, votabletype = tuple(tablesfile.iter_tables())

        self._test_datatypes_votable(votable.columns)
        self._test_datatypes_votable(votabletype.columns)

    def _test_datatypes_tap(self, cols):
        assert cols[0].datatype.content == 'BOOLEAN'
        assert cols[1].datatype.content == 'SMALLINT'
        assert cols[2].datatype.content == 'INTEGER'
        assert cols[3].datatype.content == 'BIGINT'
        assert cols[4].datatype.content == 'REAL'
        assert cols[5].datatype.content == 'DOUBLE'
        assert cols[6].datatype.content == 'TIMESTAMP'
        assert cols[7].datatype.content == 'CHAR'
        assert cols[8].datatype.content == 'VARCHAR'
        assert cols[9].datatype.content == 'BINARY'
        assert cols[10].datatype.content == 'VARBINARY'
        assert cols[11].datatype.content == 'POINT'
        assert cols[12].datatype.content == 'REGION'
        assert cols[13].datatype.content == 'CLOB'
        assert cols[14].datatype.content == 'BLOB'

    def test_datatypes_tap(self):
        tablesfile = vosi.parse_tables(
            get_pkg_data_filename("data/tables/datatypes_tap.xml"))
        tap, taptype = tuple(tablesfile.iter_tables())

        self._test_datatypes_tap(tap.columns)
        self._test_datatypes_tap(taptype.columns)

    def test_wrong_datatypes_tap(self):
        with pytest.warns(W02):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/wrong_datatypes_tap.xml"))

    def test_wrong_datatypes_votable(self):
        with pytest.warns(W02):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/wrong_datatypes_votable.xml"))

    def test_no_schemas(self):
        with pytest.warns(W14):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/no_schemas.xml"))

        with pytest.raises(W14):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/no_schemas.xml"),
                pedantic=True)

    def test_no_schema_name(self):
        with pytest.raises(E06):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/no_schema_name.xml"))

    def test_multiple_schema_names(self):
        with pytest.warns(W05):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/multiple_schema_names.xml"))

        with pytest.raises(W05):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/multiple_schema_names.xml"),
                pedantic=True)

    def test_multiple_schema_titles(self):
        with pytest.warns(W13):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_schema_titles.xml"))

        with pytest.raises(W13):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_schema_titles.xml"),
                pedantic=True)

    def test_multiple_schema_descriptions(self):
        with pytest.warns(W06):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_schema_descriptions.xml"))

        with pytest.raises(W06):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_schema_descriptions.xml"),
                pedantic=True)

    def test_multiple_schema_utypes(self):
        with pytest.warns(W09):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_schema_utypes.xml"))

        with pytest.raises(W09):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_schema_utypes.xml"),
                pedantic=True)

    def test_no_table_name(self):
        with pytest.raises(E06):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/no_table_name.xml"))

    def test_multiple_table_names(self):
        with pytest.warns(W05):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/multiple_table_names.xml"))

        with pytest.raises(W05):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/multiple_table_names.xml"),
                pedantic=True)

    def test_multiple_table_titles(self):
        with pytest.warns(W13):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_table_titles.xml"))

        with pytest.raises(W13):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_table_titles.xml"),
                pedantic=True)

    def test_multiple_table_descriptions(self):
        with pytest.warns(W06):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_table_descriptions.xml"))

        with pytest.raises(W06):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_table_descriptions.xml"),
                pedantic=True)

    def test_multiple_table_utypes(self):
        with pytest.warns(W09):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_table_utypes.xml"))

        with pytest.raises(W09):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_table_utypes.xml"),
                pedantic=True)

    def test_multiple_column_names(self):
        with pytest.warns(W05):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/multiple_column_names.xml"))

        with pytest.raises(W05):
            vosi.parse_tables(
                get_pkg_data_filename("data/tables/multiple_column_names.xml"),
                pedantic=True)

    def test_multiple_column_descriptions(self):
        with pytest.warns(W06):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_column_descriptions.xml"))

        with pytest.raises(W06):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_column_descriptions.xml"),
                pedantic=True)

    def test_multiple_column_units(self):
        with pytest.warns(W07):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_column_units.xml"))

        with pytest.raises(W07):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_column_units.xml"),
                pedantic=True)

    def test_multiple_column_ucds(self):
        with pytest.warns(W08):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_column_ucds.xml"))

        with pytest.raises(W08):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_column_ucds.xml"),
                pedantic=True)

    def test_multiple_column_utypes(self):
        with pytest.warns(W09):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_column_utypes.xml"))

        with pytest.raises(W09):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_column_utypes.xml"),
                pedantic=True)

    def test_multiple_column_datatypes(self):
        with pytest.warns(W37):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_column_datatypes.xml"))

        with pytest.raises(W37):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_column_datatypes.xml"),
                pedantic=True)

    def test_tap_size(self):
        with pytest.warns(W03):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/sizenegative.xml"))

        with pytest.raises(W03):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/sizenegative.xml"),
                pedantic=True)

    @pytest.mark.xfail
    def test_wrong_flag(self):
        with pytest.warns(W04):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/wrong_flag.xml"))

        with pytest.raises(W04):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/wrong_flag.xml"),
                pedantic=True)

    def test_multiple_fromcolumns(self):
        with pytest.warns(W10):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_fromcolumns.xml"))

        with pytest.raises(W10):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_fromcolumns.xml"),
                pedantic=True)

    def test_missing_fromcolumn(self):
        with pytest.raises(E02):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/no_fromcolumn.xml"))

    def test_multiple_targetcolumns(self):
        with pytest.warns(W11):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_targetcolumns.xml"))

        with pytest.raises(W11):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_targetcolumns.xml"),
                pedantic=True)

    def test_missing_targetcolumn(self):
        with pytest.raises(E03):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/no_targetcolumn.xml"))

    def test_multiple_targettables(self):
        with pytest.warns(W12):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_targettables.xml"))

        with pytest.raises(W12):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_targettables.xml"),
                pedantic=True)

    def test_multiple_foreignkey_descriptions(self):
        with pytest.warns(W06):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_foreignkey_descriptions.xml"))

        with pytest.raises(W06):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_foreignkey_descriptions.xml"),
                pedantic=True)

    def test_multiple_foreignkey_utypes(self):
        with pytest.warns(W09):
            vosi.parse_tables(get_pkg_data_filename(
                "data/tables/multiple_foreignkey_utypes.xml"))

        with pytest.raises(W09):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/multiple_foreignkey_utypes.xml"),
                pedantic=True)

    def test_wrong_arraysize(self):
        with pytest.raises(E01):
            vosi.parse_tables(
                get_pkg_data_filename(
                    "data/tables/wrong_arraysize.xml"))

    def test_no_table_description(self):
        """Test handling of describing tables with no description
        """
        tableset = vosi.parse_tables(
            get_pkg_data_filename(
                "data/tables/no_table_description.xml"))
        nodesc_table = tableset.get_first_table()
        assert nodesc_table.description is None

        with io.StringIO() as buf, contextlib.redirect_stdout(buf):
            nodesc_table.describe()
            output = buf.getvalue()
        assert 'No description' in output

    def test_single_table_description(self):
        """Test describing a table with a single description
        """
        tableset = vosi.parse_tables(
            get_pkg_data_filename(
                "data/tables/single_table_description.xml"))
        onedesc_table = tableset.get_first_table()
        describe_string = 'A test table with a single description'
        assert describe_string in onedesc_table.description

        with io.StringIO() as buf, contextlib.redirect_stdout(buf):
            onedesc_table.describe()
            output = buf.getvalue()
        assert describe_string in output
