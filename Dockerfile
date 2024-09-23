# Use default Python image as the base image 
FROM python:3.9-slim

# Set up working directory inside the container 
WORKDIR /crossdome-p

# Copy requirements file into container 
COPY requirements.txt . 

# Install Python3 dependencies listed in [requirements.txt]
RUN pip install --no-cache-dir -r requirements.txt 

# Copy the working directory's code into the container 
COPY . . 

# Expose ports as needed [N/A] 

# Define the default command to run main program 
# TBD 