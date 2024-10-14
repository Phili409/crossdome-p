# Use default Python image as the base image 
FROM python:3.12-slim

# Set up working directory inside the container 
WORKDIR /app

# Install system dependencies 
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    libffi-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file into container 
COPY requirements.txt . 

# Install Python3 dependencies 
RUN pip install --upgrade pip && \
    pip install -r requirements.txt 

# Copy the working directory's code into the container 
COPY . . 

# Expose ports as needed [N/A] 

# Define the default command to run main program 
# TBD 
CMD ["python", "-m", "unittest", "discover", "-s", "tests"]