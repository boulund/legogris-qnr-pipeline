import time
import sys

_logfileseparator = "----------------------------------------------------------------------\n"

class Logger:
    """
    Simple logging class with debugging capabilities. Writes both to stdout and file.
    """
    def __init__(self, logfilepath, debug=False):
        """
        Args:
            * logfilepath (str): Path to log file.

        Kwargs:
            * debug (bool, False): Whether debug statments should be written or not.
        """

        self.logfilepath = logfilepath
        self._debug = debug

    def open(self):
        """
        Open log file and start logging
        """

        t = time.asctime(time.localtime())
        try:
            self.logfile = open(self.logfilepath,'a')
            if self.logfile.tell() == 0:
                self.write("Logfile '"+self.logfilepath+"' created on "+t+"\n")
            else:
                self.write("Logging started on "+t+"\n")
        except IOError:
            print "NOTE: cannot create logfile:", arg
            print "Messages will be printed to STDOUT exclusively"

    def write(self, message):
        """
        Write message to log file and stdout.
        """
        sys.stdout.write(message)
        self.logfile.write(message)

    def writeline(self, message):
        """
        Write message to log file and stdout, followed by a newline.
        """
        print message
        self.logfile.write(message + '\n')

    def debug(self, message):
        """
        Write message iff debug was set to True at initalization.
        """
        if self._debug:
            self.write(message)

    def debugline(self, message):
        """
        Write message followed by newline iff debug was set to True at initalization.
        """
        if self._debug:
            self.writeline(message)

    def flush(self):
        """
        Flush log to file.
        """
        self.logfile.flush()

    def line(self):
        """
        Write a separator line.
        """
        self.write(_logfileseparator)

    def close(self):
        """
        Close log file.
        """
        self.logfile.close()
