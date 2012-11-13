# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

def parse_woollam_file( fname, witherror=False ):
    """ 
    Returns three lists: wavelength, n, k
    This function parses a text file of ellipsometer data from the Woollam spectroscopic ellipsometer.
    It assumes that you did a copy-pase from the fitted data into a normal text file (which I think is
    about the only way you can save the data). It also assumes that you did a search-replace on the
    data to get rid of the +- symbol before the error. The sed command to do that is:
    sed "s/$(echo -e \\0261\\c)/ /g" < messed_up_data.txt > fixed_data.txt
    where 261 is the hex code for the +- symbol (as determined by the bash command: od -c < messed_up_data.txt)
    """
    txtfile = open(fname,'rb')

    n = [] # index of refraction
    k = [] # loss coefficient
    wavelength = []
    for line in txtfile:
        if '\xb1' in line:
            # replace the +- symbol with a space
            line = line.replace( '\xb1', ' ' )
            
        if 'Thick' in line:
            pass # you forgot to comment out the thickness
        elif 'n' in line:
            # index of refraction data
            n.append( float(line.split('\t')[1].split(' ')[0]) )
            wavelength.append( float(line.split('.')[1]) )
        elif 'k' in line:
            # loss coefficient data
            k.append( float(line.split('\t')[1].split(' ')[0]) )
        
    txtfile.close()
    return wavelength, n, k

class TxtFile():
    """ Pass in a woollam-saved text file (like copy-pasted from the fit results)
    and the methods here can plot n and k. """
    def __init__( self, fname ):
        self.fname = fname
        self.wavelength, self.n, self.k = parse_woollam_file( fname )
        
    def plot_n( self, *args, **kwargs ):
        """ Plot n versus wavelength. """
        plt.plot( self.wavelength, self.n, *args, **kwargs )
        
    def plot_k( self, *args, **kwargs ):
        """ Plot k versus wavelength. """
        plt.plot( self.wavelength, self.k, *args, **kwargs )

