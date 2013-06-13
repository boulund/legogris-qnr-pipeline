import time

_logfileseparator = "----------------------------------------------------------------------"

class Logger:
    def __init__(self, logfilepath, debug=False):
        self.logfilepath = logfilepath
        self._debug = debug

    def open(self):
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
        print message
        self.logfile.write(message)

    def debug(self, message):
        if self._debug:
            self.write(message)

    def flush(self):
        self.logfile.flush()

    def line(self):
        self.write(_logfileseparator)
