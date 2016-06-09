import datetime

#----------------------------------------------------------------------------
def save_time( filename, string, addition=True ):

    now = datetime.datetime.now()
    if ( addition ): # add to an exiting file
        with open( filename, 'a') as f :
            f.write(string+now.strftime("%Y-%m-%d %H:%M:%S") + '\n')
    else : # create a new file
        with open( filename, 'w') as f :
            f.write(string+now.strftime("%Y-%m-%d %H:%M:%S") + '\n')


#----------------------------------------------------------------------------
def save_line( filename, string, addition=True ):

    if ( addition ): # add to an exiting file
        with open( filename, 'a') as f :
            f.write( string )
            f.write( '\n' )
    else : # create a new file
        with open( filename, 'w') as f :
            f.write( string )
            f.write( '\n' )
