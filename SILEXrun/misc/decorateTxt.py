import math


class decorateTxt:

    def __init__(self, le=None, pa='#'):
        self.pattern = pa
        self.length = 15
        self.adaptLength(le)

    def adaptLength(self, le=None):
        if le:
            self.length = le

    def buildPattern(self, le=None, lt=None, txt=None):
        self.adaptLength(le)
        if lt:
            length = lt
        else:
            length = self.length
        if not txt:
            txt = self.pattern
        return txt * length

    def equalPattern(self, le=None, lt=None):
        return self.buildPattern(le, lt, '=')

    def inPattern(self, le=None, lt=None):
        return self.buildPattern(le, lt, '<')

    def outPattern(self, le=None, lt=None):
        return self.buildPattern(le, lt, '>')

    def plusPattern(self, le=None, lt=None):
        return self.buildPattern(le, lt, '+')

    def minusPattern(self, le=None, lt=None):
        return self.buildPattern(le, lt, '-')

    def starPattern(self, le=None, lt=None):
        return self.buildPattern(le, lt, '*')

    def dashPattern(self, le=None, lt=None):
        return self.buildPattern(le, lt, '#')

    def adaptTxtCenter(self, fullTxt=None, le=None, lt=None):
        return self.adaptTxt(fullTxt=fullTxt,
                             le=le,
                             lt=lt,
                             pattern='>>',
                             patternB='<<')

    def adaptTxt(self,
                 fullTxt=None,
                 le=None,
                 lt=None,
                 pattern=None,
                 patternB=None):
        self.adaptLength(le)
        if lt:
            length = lt
        else:
            length = self.length
        if not pattern:
            pattern = self.pattern
        if not patternB:
            patternB = self.pattern
        # length of text
        lenTxt = len(fullTxt)
        # minimum start and end
        minStart = pattern * 2 + ' ' * 2
        minEnd = ' ' * 2 + patternB * 2
        if lenTxt + len(minStart+minEnd) < length:
            nomSpace = math.floor(
                (length - (lenTxt + len(minStart+minEnd)))/2)
            correcSpace = (length - (lenTxt + len(minStart+minEnd))) % 2
            minStart += ' ' * nomSpace
            minEnd = ' ' * nomSpace + minEnd
            if correcSpace == 1:
                minEnd = ' ' + minEnd
        return minStart + fullTxt + minEnd
