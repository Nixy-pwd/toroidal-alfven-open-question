"""Compute conservative spectral mass and duration bounds from an exact certificate."""
from fractions import Fraction as F
import argparse
from verify import verify

def consequences(error):
    frequency=F(41491,100000)
    delta=F(1,10000)
    gap=2*frequency*delta-delta**2
    mass_lower=1-(error/gap)**2
    pct_floor=(10000*mass_lower.numerator)//mass_lower.denominator
    print('Spectral weight in [0.41481,0.41501] at least',str(F(pct_floor,100)),'percent')
    def safe(periods):
        time_upper=F(44*periods,7)/frequency
        return (time_upper*error)**2 <= (2*frequency**2-error)/100
    periods=0
    while safe(periods+1):periods+=1
    print('Full ideal evolution: at most 10% relative L2 error through',periods,'whole periods')
    return periods

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('trial')
    _,err=verify(ap.parse_args().trial)
    consequences(err)
