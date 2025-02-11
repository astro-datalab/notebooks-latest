from contextlib import contextmanager

import pytest
import requests_mock


class ContextAdapter(requests_mock.Adapter):
    """
    requests_mock adapter where ``register_uri`` returns a context manager
    """
    @contextmanager
    def register_uri(self, *args, **kwargs):
        matcher = super().register_uri(*args, **kwargs)

        yield matcher

        self.remove_matcher(matcher)

    def remove_matcher(self, matcher):
        if matcher in self._matchers:
            self._matchers.remove(matcher)


@pytest.fixture(scope='function')
def mocker():
    with requests_mock.Mocker(
        adapter=ContextAdapter(case_sensitive=True)
    ) as mocker_ins:
        yield mocker_ins
