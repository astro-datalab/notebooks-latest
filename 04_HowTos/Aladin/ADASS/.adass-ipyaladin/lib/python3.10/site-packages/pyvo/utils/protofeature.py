from dataclasses import dataclass

__all__ = ['Feature']


@dataclass
class Feature:
    """
    A prototype feature implementing a standard that is currently in the process of being
    approved, but that might change as a result of the approval process. A Feature must
    have a name. Optionally, a feature may have a *url* that is displayed to the user in
    case a feature is used without the user explicitly opting in on its usage. The URL is
    expected to contain more information about the standard and its state in the approval process.
    """
    name: str
    url: str = ''
    on: bool = False

    def should_error(self):
        """
        Should accessing this feature fail?

        Returns
        -------
        bool Whether accessing this feature should result in an error.
        """
        return not self.on

    def error(self, function_name):
        """
        Format an error message when the feature is being accesses without the user having opted in its usage.

        This function will be used as a callback when an error message needs to be displayed to
        the user, with the function name that was accessed as an argument. Extensions of this
        class may have additional information to display.

        Parameters
        ----------
        function_name: str
            The name of the function associated to this feature and that the user called.

        Returns
        -------
        str: The error message to be displayed to the user.

        """
        message = (f'{function_name} is part of a prototype feature ({self.name}) that has not '
                   'been activated. For information about prototype features please refer to '
                   'https://pyvo.readthedocs.io/en/latest/utils/prototypes.html .')
        if self.url:
            message += f' For more information about the {self.name} feature please visit {self.url}.'
        message += (" To suppress this error and enable the feature use "
                    f"`pyvo.utils.activate_features('{self.name}')`")
        return message
