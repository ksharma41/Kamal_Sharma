#!/usr/bin/env python
# coding: utf-8

# In[1]:


#NY_Taxi dataset
# Reference: https://www.kaggle.com/gaborfodor/from-eda-to-the-top-lb-0-367/comments
#%matplotlib inline
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from datetime import timedelta
import datetime as dt
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [16, 10]
import seaborn as sns
#import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans


# In[2]:


#read data once
train_df = pd.read_csv('train.csv')
test_df = pd.read_csv('test.csv')


# In[3]:


#dataframe.shape[0] or dataframe.shape[1] gives the length along the y and x direction respectively.
print ("train_df has %s rows and %s coloumns" % (train_df.shape[0],train_df.shape[1]))
print ("train_df has %s rows and %s coloumns" % (test_df.shape[0],test_df.shape[1]))
train_df.head()


# In[4]:


#lets covert the date and time in train and test files into datetime stamp
train_df['pickup_datetime'] = pd.to_datetime(train_df.pickup_datetime)
train_df['dropoff_datetime'] = pd.to_datetime(train_df.dropoff_datetime)

test_df['pickup_datetime'] = pd.to_datetime(test_df.pickup_datetime)

#separating pickup date and pickup time
train_df.loc[:,'pickup_date'] = train_df['pickup_datetime'].dt.date
test_df.loc[:,'pickup_date'] = test_df['pickup_datetime'].dt.date

train_df.loc[:,'pickup_time'] = train_df['pickup_datetime'].dt.time
test_df.loc[:,'pickup_time'] = test_df['pickup_datetime'].dt.time

#Converting 'Y' in 'store_and_forward' column to Boolean 1 and 0.
train_df['store_and_fwd_flag'] = 1 * (train_df.store_and_fwd_flag.values == 'Y')
test_df['store_and_fwd_flag'] = 1 * (test_df.store_and_fwd_flag.values == 'Y')

#output shows that the string values in pickup_datetime have been converted to timestamps.


# In[5]:


#Check whether trip durations in the original training data are correct or not. we calculate
# difference between given duration with the one we will calculate.
train_df['calculated_duration'] = (train_df['dropoff_datetime'] - train_df['pickup_datetime']).map(lambda x: x.total_seconds())
#Check !!
if (np.array_equal(train_df['calculated_duration'], train_df['trip_duration'])): 
    print ('All Iz Well !!') 
else: print("Durations not equal !!")


# In[6]:


# Now that the data is all good. Lets have a look at the distribution of trip_durations.
#Histogram
train_df['log_trip_duration'] = np.log(train_df['trip_duration'])
plt.hist((train_df['log_trip_duration'].values) , bins = 'auto')
plt.title('Histogram of logarithm of trip_duration')
plt.ylabel('Frequency of instances')
plt.xlabel('log(trip_duration)')
plt.show()


# In[7]:


#Let's see if there are any partially missing values in test and train data 
print ("All columns have the same number of non 'Nan' values. Therefore, no missing values. ")
print train_df.count()
print test_df.count()


# The number of elements in each column is the same for both train and test data. This means there are no 'Nan' values (missing values).

# #  Check for time frame overlap between train and test data
# We now group the trips based on the pickup dates and plot the number of trips for each day.
# It is important to see whether the train data and test data  belong to the same time frames and geographical region.
# e.g. We cannot find trustable models detailed anough to predict traffic in SanFrancisco based on New York data.
# Or the data from March to create good models for traffic in December.

# In[8]:


#Check for time frame Overlap
plt.plot(train_df.groupby('pickup_date').count()[['id']],'.-',color = 'red', label = 'Train Data')
plt.plot(test_df.groupby('pickup_date').count()[['id']],'.-',color = 'blue', label = 'Train Data')
plt.title('Trips on different calendar days')
plt.ylabel('Number of trips')
plt.xlabel('Calendar Dates')
plt.legend(loc=0)
print ("The time frame for the data collected overlaps.")


# # Check for geographical region overlap - Visualization
# For the same reason as above, we check the geographical overlap between train and test data. The taxis in New York that are in the same neighborhood wil have to go through the same traffic and same city blocks. Let's see how our data is distributed. We will now plot the pickup and dropoff positions treating latitudes and longitudes as $x$ and $y$ coordinates.

# In[9]:


# check for geographical overlap
west_border_long = min(min(train_df['pickup_longitude'].values),min(test_df['pickup_longitude'].values))
east_border_long = max(max(train_df['pickup_longitude'].values),max(test_df['pickup_longitude'].values))
north_border_lat = max(max(train_df['pickup_latitude'].values),max(test_df['pickup_latitude'].values))
south_border_lat = min(min(train_df['pickup_latitude'].values),min(test_df['pickup_latitude'].values))

print ("The geographical edges of the available training data are latitude %f N to %f N " % (south_border_lat, north_border_lat)),
print ("and longitude %f E to %f E " % (west_border_long, east_border_long))


# In[26]:


plt.hist((test_df['pickup_latitude'].values[:len(test_df['pickup_latitude'].values):100]) , bins = 'auto')
plt.title('Distribution of Pick-up Latitudes ')
plt.xlabel('Latitude')
plt.ylabel('Number of Trips')
plt.show()

plt.hist((test_df['pickup_longitude'].values[:len(test_df['pickup_longitude'].values):100]) , bins = 'auto')
plt.title('Distribution of Pick-up Longitudes ')
plt.xlabel('Longitude')
plt.ylabel('Number of Trips')
plt.show()


# In[43]:


#city_long_border = (west_border_long, east_border_long)
#city_lat_border = (south_border_lat,north_border_lat)

m = 10;
n = 1 * m;
city_long_border = (-74.05, -73.8)
city_lat_border = (40.6, 40.85)
fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True)
ax[0].scatter(train_df['pickup_longitude'].values[:1458644:n], train_df['pickup_latitude'].values[:1458644:n],
              color='blue', s=1, label='train', alpha=0.015)
ax[1].scatter(test_df['pickup_longitude'].values[:625134:m], test_df['pickup_latitude'].values[:625134:m],
              color='green', s=1, label='test', alpha=0.015)
fig.suptitle('Train and test area complete overlap.')
ax[0].legend(loc=1)
ax[0].set_ylabel('latitude')
ax[0].set_xlabel('longitude')
ax[1].set_xlabel('longitude')
ax[1].legend(loc=0)

plt.ylim(city_lat_border)
plt.xlim(city_long_border)
plt.show()


# In[ ]:




