from astroquery.irsa import Irsa
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy import units as u

def setup(ra_min=23.43, ra_max=23.44, dec_min=30.71, dec_max=30.72):
        ###plots unWISE sources
    
    Irsa.ROW_LIMIT = 100000
    #table = Irsa.query_region("m33", catalog="unwise_2019", spatial="Polygon", radius=0.25*u.deg)
    table = Irsa.query_region("m33", catalog="unwise_2019", spatial="Polygon",
            polygon=[SkyCoord(ra=ra_min,dec=dec_min,unit=(u.deg,u.deg),frame='icrs'),
                     SkyCoord(ra=ra_max,dec=dec_min,unit=(u.deg,u.deg),frame='icrs'),
                     SkyCoord(ra=ra_max,dec=dec_max,unit=(u.deg,u.deg),frame='icrs'),
                     SkyCoord(ra=ra_min,dec=dec_max,unit=(u.deg,u.deg),frame='icrs')])

    fig = plt.figure()
    fig.set_figwidth(8)
    fig.set_figheight(6.5)
    ax = plt.subplot(111)
    ax.set_xlabel('RA', fontsize=20)
    ax.set_ylabel('Dec', fontsize=20)

    ax.scatter(table['ra'],table['dec'],s=10,c='black')

    ax.set_xlim(ax.get_xlim()[::-1])
    
    x,y = table['ra'].round(4),table['dec'].round(4)
    
    for i_x, i_y in zip(x,y):
        plt.text(i_x, i_y, '({}, {})'.format(i_x, i_y))  

        
def sources():
      interact(setup, ra_min=(23.43, 23.6, 0.01), ra_max=(23.431, 23.6, 0.01), 
         dec_min=(30.71, 30.81, 0.01), dec_max=(30.711, 30.81, 0.01))