class Colorcodes(object):
    """
    Provides ANSI terminal color codes which are gathered via the ``tput``
    utility. That way, they are portable. If there occurs any error with
    ``tput``, all codes are initialized as an empty string.

    The provides fields are listed below.

    Control:

    - bold
    - reset

    Colors:

    - blue
    - green
    - orange
    - red

    :license: MIT
    """
    def __init__(self):
        import subprocess

        try:
            self.bold = subprocess.check_output("tput bold".split())
            self.reset = subprocess.check_output("tput sgr0".split())

            self.blue = subprocess.check_output("tput setaf 4".split())
            self.green = subprocess.check_output("tput setaf 2".split())
            self.orange = subprocess.check_output("tput setaf 3".split())
            self.red = subprocess.check_output("tput setaf 1".split())
        except subprocess.CalledProcessError as e:
            self.bold = b""
            self.reset = b""

            self.blue = b""
            self.green = b""
            self.orange = b""
            self.red = b""

