# CD/Dockerfile
FROM ubuntu:20.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    python3 \
    python3-pip \
    libboost-all-dev \
    libeigen3-dev

# Set working directory
WORKDIR /app

# Copy C++ source code
COPY src/calculators /app/src/calculators
COPY src/utils /app/src/utils
COPY CMakeLists.txt /app/

# Build C++ library
RUN mkdir -p build && cd build && \
    cmake .. && \
    make

# Copy Python API code
COPY src/api /app/src/api
COPY requirements.txt /app/

# Install Python dependencies
RUN pip3 install -r requirements.txt

# Set environment variables
ENV PYTHONPATH=/app

# Expose port
EXPOSE 5000

# Run API server
CMD ["python3", "src/api/app.py"]