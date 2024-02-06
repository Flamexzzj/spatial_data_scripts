# -*- coding: utf-8 -*-


import os
from datetime import datetime
from glob import glob

import numpy as np
import rasterio
from shapely.geometry import Polygon, mapping
from tqdm import tqdm

import pystac
from pystac.extensions.eo import Band, EOExtension
from pystac.extensions.scientific import ScientificExtension
from pystac.extensions.scientific import Publication

from pystac.provider import Provider
from pystac.provider import ProviderRole

import pystac.provider 

from pystac import Provider
import geopandas as gpd
from shapely.geometry import shape
from shapely.geometry import box
import pandas as pd


def get_bbox_and_footprint(raster_uri):
	with rasterio.open(raster_uri) as ds:
		bounds = ds.bounds
		bbox = [bounds.left, bounds.bottom, bounds.right, bounds.top]
		footprint = Polygon([
			[bounds.left, bounds.bottom],
			[bounds.left, bounds.top],
			[bounds.right, bounds.top],
			[bounds.right, bounds.bottom]
		])
		
		return (bbox, mapping(footprint))
	
def get_bbox_and_footprint_json(geojson_uri):
    # Load the GeoJSON file
    gdf = gpd.read_file(geojson_uri)

    # Assuming the GeoJSON contains at least one feature
    # and you are interested in the first feature
    feature = gdf.iloc[0].geometry

    # Calculate the bounding box
    bounds = feature.bounds
    bbox = [bounds[0], bounds[1], bounds[2], bounds[3]]

    # Create a footprint from the bounding box
    footprint = box(*bounds)

    return (bbox, mapping(footprint))

events = ['US-Texas', 'US-Carolina', 'Ghana', 'Nepal', 'US-Dakota', 'Spain', 'US-Nebraska', 'Uzbekistan', 
		  'Nigeria', 'US-Alabama', 'Somalia', 'Bolivia', 'Colombia', 'US-Oklahoma', 'Bangladesh', 
		  'US-Kansas', 'US-Arkansas', 'Cambodia', 'Paraguay']
folder_path = "/Users/zhijiezhang/Current_Projects/CSDA_Project/Data_done/Zhijie_FloodPlanet_2023"
sensors = ['L8','labels', 'PS', 'S1', 'S2']
img_test = "/Users/zhijiezhang/Current_Projects/CSDA_Project/Data_done/Zhijie_FloodPlanet_2023/Bolivia/PS/BOL_1311.tif"
bound, _ = get_bbox_and_footprint(img_test)

license = 'CC-BY-SA-4.0'
items_mission_description = 'This FloodPlanet dataset is a manually labeled inundation dataset that contains 19 flood events. It contains manual inundation labels based on PlanetScope, modified 8-bit version of PlanetScope data, corresponding Sentinel-1&2 imagery and Landsat-8 imagery.'
#Provider information
provider1 = Provider(name='Zhijie(JJ) Zhang',
                       description="Postdoctoral Research Associate, University of Arizona, Social [Pixel] Lab",
                       roles=[ProviderRole.PRODUCER, ProviderRole.PROCESSOR],
                       url='https://scholar.google.com/citations?hl=zh-CN&user=1Uyru5oAAAAJ')
provider2 = Provider(name='Alex Melancon',
                       description="Research Associate, University of Alabama Huntsville",
                       roles=[ProviderRole.PROCESSOR],
                       url='https://www.linkedin.com/in/alexander-melancon-865165154')
provider3 = Provider(name='Beth Tellman',
                       description="Assistant Professor, University of Arizona, Social [Pixel] Lab",
                       roles=[ProviderRole.PRODUCER],
                       url='https://beth-tellman.github.io/index.html')
provider4 = Provider(name='Jonathan Giezendanner',
                       description="Postdoctoral Researcher University of Arizona, Social [Pixel] Lab",
                       roles=[ProviderRole.PROCESSOR],
                       url='https://scholar.google.ch/citations?hl=en&user=cNxEnK4AAAAJ&view_op=list_works&sortby=pubdate')
provider5 = Provider(name='Rohit Mukherjee',
                       description="Postdoctoral Researcher University of Arizona, Social [Pixel] Lab",
                       roles=[ProviderRole.PROCESSOR],
                       url='https://rohitmukherjee.space/')
provider1_dict = provider1.to_dict()
provider2_dict = provider2.to_dict()
provider3_dict = provider3.to_dict()
provider4_dict = provider4.to_dict()
provider5_dict = provider5.to_dict()

contact_email = 'zhijiezhang@arizona.edu'

catalog = pystac.Catalog(id='FloodPlanet Manually Labeled Inundation Dataset', description='STAC catalog for FloodPlanet dataset', title='FloodPlanet Manually Labeled Inundation Dataset')

catalog.extra_fields['license'] = 'CC-BY-SA-4.0'
catalog.extra_fields['Temporal extent'] = '2017-03-14 ~ 2020-05-10'
catalog.extra_fields['Contact'] = contact_email
catalog.extra_fields['Providers'] = [provider1_dict, provider2_dict, provider3_dict, provider4_dict, provider5_dict]
catalog.extra_fields['Mission'] = items_mission_description
# extent_coords = [[-180, -90, 180, 90]]
spatial_extent = pystac.SpatialExtent(bboxes=[bound])
temporal_all = ['2017-03-14T00:00:00Z','2020-05-10T00:00:00Z']
dates_all = [datetime.strptime(d, '%Y-%m-%dT%H:%M:%SZ') for d in temporal_all]
temporal_extent_all = pystac.TemporalExtent(intervals=[dates_all])

Planet_bands = [Band.create(name='Band 1', description='Blue'),
				Band.create(name='Band 2', description='Green'),
				Band.create(name='Band 3', description='Red'),
				Band.create(name='Band 4', description='NIR')
				]
S2_bands = [Band.create(name='Band 1', description='blue'),
				Band.create(name='Band 2', description='green'),
				Band.create(name='Band 3', description='red'),
				Band.create(name='Band 4', description='Red Edge 1'),
				Band.create(name='Band 5', description='Red Edge 2'),
				Band.create(name='Band 6', description='Red Edge 3'),
				Band.create(name='Band 7', description='NIR'),
				Band.create(name='Band 8', description='Red Edge 4'),
				Band.create(name='Band 9', description='SWIR 1'),
				Band.create(name='Band 10', description='SWIR 2'),
				]
L8_bands = [Band.create(name='Band 1', description='Ultra Blue'),
				Band.create(name='Band 2', description='Blue'),
				Band.create(name='Band 3', description='Green'),
				Band.create(name='Band 4', description='Red'),
				Band.create(name='Band 5', description='Near Infrared'),
				Band.create(name='Band 6', description='Shortwave Infrared 1'),
				Band.create(name='Band 7', description='Shortwave Infrared 2'),
				]
S1_bands = [Band.create(name='Band 1', description='VV'),
			Band.create(name='Band 2', description='VH')
			]
label_bands = [Band.create(name='Value 0', description='No Data'),
			Band.create(name='Value 1', description='Not Water'),
			Band.create(name='Value 2', description='Inundation')
			]

date_file_path = '/Users/zhijiezhang/Desktop/csda-sensor-dates.xlsx'
df = pd.read_excel(date_file_path)



# Extract the date information for the specified event and sensor
date_info1 = df[df['Events'] == event]['PS'].iloc[0]
date_info2 = df[df['Events'] == event]['S2'].iloc[0]

# Define the temporal extent (start and end time)
temporal_extent_event = pystac.TemporalExtent(
    intervals=[[date_info1, date_info2]]
)

for event in tqdm(events):
	event_json = '/Users/zhijiezhang/Current_Projects/CSDA_Project/Data_done/Event_json'+ '/' + event + '.geojson'
	event_bound, _ = get_bbox_and_footprint_json(event_json)
	spatial_extent_event = pystac.SpatialExtent(bboxes=[event_bound])##
	date_info1 = pd.to_datetime(df[df['Events'] == event]['PS'].iloc[0])
	date_info2 = pd.to_datetime(df[df['Events'] == event]['S2'].iloc[0])

	# Define the temporal extent (start and end time)
	temporal_extent_event = pystac.TemporalExtent(
		intervals=[[date_info1, date_info2]]
	)
	collection_extent = pystac.Extent(spatial=spatial_extent_event, temporal=temporal_extent_event)

	event_first_collection = pystac.Collection(id=event, title='Flood event occured in '+event, description='Folders under this event', extent=collection_extent)
	event_first_collection.license = license
	event_first_collection.extra_fields['Contact'] = contact_email
	event_first_collection.extra_fields['Providers'] = [provider1_dict, provider2_dict, provider3_dict, provider4_dict]
	event_first_collection.extra_fields['Mission'] = items_mission_description

	event_first_collection.catalog = catalog
	bboxes = []
	for sensor in sensors:
		# date_info3 = pd.to_datetime(df[df['Events'] == event][sensor].iloc[0], format='%m-%d-%Y')
		date_info3 = pd.to_datetime(df[df['Events'] == event][sensor].iloc[0], errors='coerce', format='%m-%d-%Y')
		# Define the temporal extent (start and end time)
		temporal_extent_event_sensor = pystac.TemporalExtent(
			intervals=[[date_info3, date_info3]]
		)
		second_collection_extent = pystac.Extent(spatial=spatial_extent_event, temporal=temporal_extent_event_sensor)
		sensor_second_collection = pystac.Collection(id=sensor, title=sensor + ' imagery', description = sensor + ' imagery', extent=second_collection_extent)
		sensor_second_collection.license = license
		sensor_second_collection.extra_fields['Contact'] = contact_email
		# sensor_second_collection.extra_fields['Providers'] = [provider_dict]
		sensor_second_collection.extra_fields['Mission'] = items_mission_description
		sensor_second_collection.catalog = catalog
		sensor_second_collection.parent = event_first_collection
		img_path = os.path.join(folder_path, event, sensor)
		img_list = glob(img_path + "/*.tif")
		for img in img_list:
			print(img)
			bbox, footprint = get_bbox_and_footprint(img)
			
			
			item_date_string = df[df['Events'] == event]['PS'].iloc[0]
			extracted_item_datetime = pd.to_datetime(item_date_string, format='%m-%d-%Y')
			item_id = img.split('/')[-1]
			# print(item_id)
			item = pystac.Item(id=item_id.split('.')[0], collection = sensor_second_collection, geometry=footprint, bbox=bbox, datetime=extracted_item_datetime, properties={})

			if sensor == 'PS':
				
				bboxes.append(bbox)

				item.properties['Event'] = event
				item.properties['Chip_ID'] = item_id.split('.')[0]
				item.properties['Image source'] = 'PlanetScope'
				item.properties['license'] = 'CC-BY-SA-4.0'
				item.properties['Contact'] = contact_email
				# item.properties['Providers'] = [provider_dict]
				item.properties['Mission'] = items_mission_description

				# item.properties['Data Description'] = 'Band 1: Red band. Band 2: Green band. Band 3: Blue band. Band 4: NIR band.'
				item.add_asset(key=sensor, asset = pystac.Asset(href=img, media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='labels', asset = pystac.Asset(href=img.replace('PS','labels'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Sentinel-1 data', asset = pystac.Asset(href=img.replace('PS','S1'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Sentinel-2 data', asset = pystac.Asset(href=img.replace('PS','S2'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Landsat-8 data', asset = pystac.Asset(href=img.replace('PS','L8'), media_type=pystac.MediaType.GEOTIFF))
				item.common_metadata.gsd = 3
				eo = EOExtension.ext(item, add_if_missing=True)
				eo.apply(bands=Planet_bands)
				
			if sensor == 'S1':
				
				bboxes.append(bbox)

				item.properties['Event'] = event
				item.properties['Chip_ID'] = item_id.split('.')[0]
				item.properties['Image source'] = 'Sentinel-1'
				item.properties['license'] = 'CC-BY-SA-4.0'
				item.properties['Contact'] = contact_email
				# item.properties['Providers'] = [provider_dict]
				item.properties['Mission'] = items_mission_description
				# item.properties['Data Description'] = 'Band 1: Red band. Band 2: Green band. Band 3: Blue band. Band 4: NIR band.'
				item.add_asset(key=sensor, asset = pystac.Asset(href=img, media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='labels', asset = pystac.Asset(href=img.replace('S1','labels'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding PlanetScope data', asset = pystac.Asset(href=img.replace('S1','PS'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Sentinel-2 data', asset = pystac.Asset(href=img.replace('S1','S2'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Landsat-8 data', asset = pystac.Asset(href=img.replace('S1','L8'), media_type=pystac.MediaType.GEOTIFF))
				item.common_metadata.gsd = 10
				eo = EOExtension.ext(item, add_if_missing=True)
				eo.apply(bands=S1_bands)

			if sensor == 'S2':
				
				bboxes.append(bbox)

				item.properties['Event'] = event
				item.properties['Chip_ID'] = item_id.split('.')[0]
				item.properties['Image source'] = 'Sentinel-2'
				item.properties['license'] = 'CC-BY-SA-4.0'
				item.properties['Contact'] = contact_email
				# item.properties['Providers'] = [provider_dict]
				item.properties['Mission'] = items_mission_description
				# item.properties['Data Description'] = 'Band 1: Red band. Band 2: Green band. Band 3: Blue band. Band 4: NIR band.'
				item.add_asset(key=sensor, asset = pystac.Asset(href=img, media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='labels', asset = pystac.Asset(href=img.replace('S2','labels'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Sentinel-1 data', asset = pystac.Asset(href=img.replace('S2','S1'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding PlanetScope data', asset = pystac.Asset(href=img.replace('S2','PS'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Landsat-8 data', asset = pystac.Asset(href=img.replace('S2','L8'), media_type=pystac.MediaType.GEOTIFF))
				item.common_metadata.gsd = 10
				eo = EOExtension.ext(item, add_if_missing=True)
				eo.apply(bands=S2_bands)

			if sensor == 'L8':
				
				bboxes.append(bbox)

				item.properties['Event'] = event
				item.properties['Chip_ID'] = item_id.split('.')[0]
				item.properties['Image source'] = 'Landsat-8'
				item.properties['license'] = 'CC-BY-SA-4.0'
				item.properties['Contact'] = contact_email
				# item.properties['Providers'] = [provider_dict]
				item.properties['Mission'] = items_mission_description
				# item.properties['Data Description'] = 'Band 1: Red band. Band 2: Green band. Band 3: Blue band. Band 4: NIR band.'
				item.add_asset(key=sensor, asset = pystac.Asset(href=img, media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='labels', asset = pystac.Asset(href=img.replace('L8','labels'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Sentinel-1 data', asset = pystac.Asset(href=img.replace('L8','S1'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding PlanetScope data', asset = pystac.Asset(href=img.replace('L8','PS'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Sentinel-1 data', asset = pystac.Asset(href=img.replace('L8','S2'), media_type=pystac.MediaType.GEOTIFF))
				item.common_metadata.gsd = 30
				eo = EOExtension.ext(item, add_if_missing=True)
				eo.apply(bands=L8_bands)

			if sensor == 'labels':
				bboxes.append(bbox)

				item.properties['Event'] = event
				item.properties['Chip_ID'] = item_id.split('.')[0]
				item.properties['Image source'] = 'Hand label'
				item.properties['Data Description'] = '0: NoData. 1: Not Water. 2: Low confidence Water. 3:Water'
				item.properties['license'] = 'CC-BY-SA-4.0'
				item.properties['Contact'] = contact_email
				# item.properties['Providers'] = [provider_dict]
				item.properties['Mission'] = items_mission_description
				item.add_asset(key='images', asset = pystac.Asset(href=img.replace('labels','L8'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Sentinel-1 data', asset = pystac.Asset(href=img.replace('labels','S1'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding PlanetScope data', asset = pystac.Asset(href=img.replace('labels','PS'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key='Corresponding Sentinel-1 data', asset = pystac.Asset(href=img.replace('labels','S2'), media_type=pystac.MediaType.GEOTIFF))
				item.add_asset(key=sensor, asset = pystac.Asset(href=img, media_type=pystac.MediaType.GEOTIFF))
				item.common_metadata.gsd = 3
				eo = EOExtension.ext(item, add_if_missing=True)
				eo.apply(bands=label_bands)
			sensor_second_collection.add_item(item)

		event_first_collection.add_child(sensor_second_collection)

	event_extent = pystac.SpatialExtent(bboxes=[np.array(bboxes)[:,0].min(),np.array(bboxes)[:,1].min(),np.array(bboxes)[:,2].max(),np.array(bboxes)[:,3].max()])
	event_first_collection.extent.spatial = event_extent

	catalog.add_child(event_first_collection)
	# catalog.add_child(sensor_second_collection)

	


# Set self_href for each Collection
for collection in catalog.get_all_collections():
	print(collection.get_links('self'))
	for i in range(len(collection.get_links('self'))):
		collection.set_self_href(collection.get_links('self')[i].get_href())
		


catalog.normalize_hrefs(os.path.join(folder_path, 'stac_catalog'))
catalog.make_all_asset_hrefs_relative()
catalog.save(catalog_type=pystac.CatalogType.SELF_CONTAINED)
