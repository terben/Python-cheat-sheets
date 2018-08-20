# This file contains the essential commands from the
# astropquery.vizier cheat sheet notebook as a Python script

# necessary modules for this sheet:
import astropy.units as au
import astropy.coordinates as ac
import astroquery.vizier as av

# obtain a list of available Vizier catalogues with a reference to Gaia DR2
cat_list = av.Vizier.find_catalogs("Gaia DR2")
print({k:v.description for k,v in cat_list.items()})

# define a sky-coordinate of our position and make a first query:
sky_coord = ac.SkyCoord(ra=210.0, dec=-10.0, unit=(au.deg, au.deg))

# do a first query with a radius of 0.05 deg around our position
# with 'all' quantities available.
#
# The '**' for the columns means 'all available quantities'
gaia_query = av.Vizier(columns = ["**"])

first_query = gaia_query.query_region(sky_coord, radius=0.05*au.deg,
                                      catalog='I/345')

#first_query = gaia_query.query_region(sky_coord,
#                                      width=0.05 * au.deg,
#                                      height=0.05 * au.deg,
#                                      catalog='I/345')
# print(first_query)

print(first_query['I/345/gaia2'])

# Note: Please run the notebook yourself to see the output of the
# following cell. It would be too big and not really useful here.
print(first_query['I/345/gaia2'].info)

# define query for the columns we want:
gaia_query = av.Vizier(columns = ["RAJ2000", "DEJ2000",
                                  "e_RAJ2000", "e_DEJ2000", "Gmag", "e_Gmag"])

# the following line is necessary that we get *all* sources and not only
# the first 50:
gaia_query.ROW_LIMIT = -1
query_result = gaia_query.query_region(sky_coord, radius=0.5 * au.deg,
                                       catalog='I/345')

print(query_result['I/345/gaia2'])

# convert error quantities to units of degrees (are in milliarcseconds)
query_result['I/345/gaia2']['e_RAJ2000'] = \
  query_result['I/345/gaia2']['e_RAJ2000'].to(au.deg)
query_result['I/345/gaia2']['e_DEJ2000'] = \
  query_result['I/345/gaia2']['e_DEJ2000'].to(au.deg)

print(query_result['I/345/gaia2'])

# writes an ASCII catalogue to your disk
query_result['I/345/gaia2'].write('gaia.txt',
                                  format='ascii.commented_header',
                                  overwrite=True)
# writes an 'extended csv-catalogue' to your disk. This ASCII-format also
# contains the catalogues meta-data (units and column description) in
# a standardized form.
query_result['I/345/gaia2'].write('gaia.csv',
                                  format='ascii.ecsv',
                                  overwrite=True)
