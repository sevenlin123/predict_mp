import pandas as pd
import numpy as np
import sys, os
#from propagate_lite import propagate
from propagate import propagate_lite as propagate
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.time import Time
import time as TT


packed_date = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9,'A':10,
               'B':11, 'C':12, 'D':13, 'E':14, 'F':15, 'G':16, 'H':17, 'I':18, 'J':19,
               'K':20, 'L':21, 'M':22, 'N':23, 'O':24, 'P':25, 'Q':26, 'R':27, 'S':28,
               'T':29, 'U':30, 'V':31} 

class known_object:
    def __init__(self):
        #self.knowns = open('Distant.txt').readlines()
        self.knowns = open('MPCORB.DAT').readlines()
        self.name = None
        self.a = None
        self.e = None
        self.i = None
        self.w = None
        self.W = None
        self.M = None
        self.epoch = None
        self.get_orb()
    
    def get_orb(self):
        name, M, w, W, i, e, a, H, epoch = [],[],[],[],[],[],[],[],[]
        for o in self.knowns:
            try:
                H.append(float(o[8:12]))
                M.append(float(o[26:35]))
                w.append(float(o[37:46]))
                W.append(float(o[48:57]))
                i.append(float(o[59:68]))
                e.append(float(o[70:79]))
                a.append(float(o[92:103]))
                century = packed_date[o[20]]*100
                year = century + int(o[21:23])
                month = packed_date[o[23]]
                day = packed_date[o[24]]
                name.append(o[0:7])
                time = ['{}-{}-{}'.format(year, month, day)]
                t = Time(time, format='isot', scale='utc')
                epoch.append(t.jd[0])
            except ValueError:
                pass
        
        #print(a[0],e[0],i[0],w[0],W[0],M[0],epoch[0],name[0])
        self.name = np.array(name)
        self.a = np.array(a)
        self.e = np.array(e)
        self.i = np.array(i)*np.pi/180.
        self.w = np.array(w)*np.pi/180.
        self.W = np.array(W)*np.pi/180.
        self.M = np.array(M)*np.pi/180.
        self.H = np.array(H)
        self.epoch = np.array(epoch)
        
    def output_csv(self):
        d = {'name':self.name, 'a':self.a, 'e':self.e, 'i':self.i, 'w':self.w, 'W':self.W,
                'M':self.M, 'H':self.H, 'epoch':self.epoch}
        df = pd.DataFrame(data=d)
        df.to_csv('MPCORB.csv', index=False)

def gen_csv():
    known = known_object()
    known.output_csv()

def predict(pointing):
    field, ra, dec, mjd = pointing.split()
    jd = float(mjd) + 2400000.5
    ra = float(ra)
    dec = float(dec)
    p=propagate(np.array(known.a), np.array(known.e), np.array(known.i), np.array(known.w), np.array(known.W), np.array(known.M), np.array(known.epoch), np.zeros(len(known.a))+jd, helio=True)
    mag = known.H + 5*np.log10(p.r*(p.delta))
    ra_matched = abs(p.ra*180/np.pi - ra) < 0.17
    dec_matched = abs(p.dec*180/np.pi - dec) < 0.1
    matched = ra_matched*dec_matched
    if matched.sum() != 0:
        print(field, jd, list(known.name[matched]), p.ra[matched]*180/np.pi, p.dec[matched]*180/np.pi, list(mag[matched]))
    
def main():
    #gen_csv()
    global known
    known = pd.read_csv('MPCORB.csv')
    pointings = open(sys.argv[1]).readlines()
    predict(pointings[0])
    #list(map(predict, pointings))

if __name__ == '__main__':
    main()		  
		
