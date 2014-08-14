import cPickle
import kasey_utils as kc

class Scan():
    def __init__( self, fname=None ):
        if fname is None:
            # we're building a raster object as the data comes in
            return self
        
        if not fname.endswith('.RPK'):
            raise ValueError( 'Can only create a Scan class object from a raster pickle file; fname=%s' % fname )

        self.fname = fname
        self.read_rpk( fname )
        self.ax = None
    

    def colormap( self ):
        """ plot a colormap of the full raster scan
        """
        x = kc.centers_to_corners( self.xpoints )
        y = kc.centers_to_corners( self.ypoints )


    def read_rpk( self, fname ):
        """ read in the raster pickle
        """
        with fopen( fname ) as fobj:
            data = cPickle.load( fobj )
        
        self.data = data
        self.xpoints = data['xpoints']
        self.ypoints = data['ypoints']


    def set_axes( self, axes ):
        """ this allows you to plot the data to a particular axes
        """
        self.ax = axes
 