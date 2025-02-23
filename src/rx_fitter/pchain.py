'''
Module holding PChain class
'''
#----------------------------------
class PChain:
    '''
    Class meant to represent a decay chain
    '''
    #----------------------------------
    def __init__(self, pid, mid, gmid, ggmid):
        self._TID          = pid
        self._MOTHER_TID   = mid
        self._GMOTHER_TID  = gmid
        self._GGMOTHER_TID = ggmid
    #----------------------------------
    def ID(self):
        return self._TID

    def MID(self):
        return self._MOTHER_TID

    def GMID(self):
        return self._GMOTHER_TID

    def GGMID(self):
        return self._GGMOTHER_TID
    #----------------------------------------------------
    def MatchDecay(self, l_dec_id):
        if len(l_dec_id) == 1:
            return self._TID ==      l_dec_id[0]
        if len(l_dec_id) == 2:
            return  self._TID == abs(l_dec_id[0]) and self._MOTHER_TID == abs(l_dec_id[1])
        if len(l_dec_id) == 3:
            return  self._TID == abs(l_dec_id[0]) and self._MOTHER_TID == abs(l_dec_id[1]) and self._GMOTHER_TID == abs(l_dec_id[2])
        if len(l_dec_id) == 4:
            return  self._TID == abs(l_dec_id[0]) and self._MOTHER_TID == abs(l_dec_id[1]) and self._GMOTHER_TID == abs(l_dec_id[2]) and self._GGMOTHER_TID == abs(l_dec_id[3]);

        return False
    #----------------------------------------------------
    def HasInChain(self, ID):
        _return = False

        if( self._MOTHER_TID   == ID):
            _return = True
        if( self._GMOTHER_TID  == ID):
            _return = True
        if( self._GGMOTHER_TID == ID):
            _return = True

        return _return
    #----------------------------------------------------
    def MatchUpstream(self, IDFirstDau, HeadPart):
        _return = False
        if(  self._MOTHER_TID   == IDFirstDau and self._GMOTHER_TID  == HeadPart):
            _return = True
        if(  self._GMOTHER_TID  == IDFirstDau and self._GGMOTHER_TID == HeadPart):
            _return = True
        if(  self._GGMOTHER_TID == IDFirstDau):
            _return = True

        return _return
    #----------------------------------------------------
    def MatchID(self, iD):
        return self._TID == abs(iD)

    def MatchMother(self, iD):
        return self._MOTHER_TID == abs(iD)

    def MatchGMother(self, iD):
        return self._GMOTHER_TID == abs(iD)

    def MatchGGMother(self, iD):
        return self._GGMOTHER_TID == abs(iD)
#----------------------------------------------------
