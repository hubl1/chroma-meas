#######################################################################
# 📦  builder  ——  基于现成 chroma-cpu 镜像来编译你的插件
#######################################################################
FROM ghcr.io/hubl1/chroma-cpu:latest AS builder
# ← 你的基础镜像 tag

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build git \
        libopenmpi-dev openmpi-bin\ 
        libhdf5-dev libpng-dev zlib1g-dev libxml2-dev&& \
    rm -rf /var/lib/apt/lists/*

ENV CC=mpicc
ENV CXX=mpicxx

WORKDIR /src

# ---- 拷贝源码 ----
COPY source/ ./source
COPY source/CMakeLists.txt ./CMakeLists.txt

# ---- 编译 ----
RUN cmake -B build \
          -DCMAKE_BUILD_TYPE=Release \
          -DChroma_DIR=/opt/chroma/lib/cmake/Chroma \
          -DCMAKE_PREFIX_PATH=/opt/chroma && \
    cmake --build build -j$(nproc) && \
    cmake --install build --prefix /opt/chroma-meas

#######################################################################
# 🚀  runtime  ——  最终镜像：基础 + 插件
#######################################################################
FROM ghcr.io/hubl1/chroma-cpu:latest
WORKDIR /meas
#COPY source/ ./source
COPY *.xml .
COPY *.sh .
COPY *.md .

# 拷贝插件安装树
COPY --from=builder /opt/chroma-meas /opt/chroma-meas

# 更新 PATH / LD_LIBRARY_PATH
ENV PATH=/opt/chroma-meas/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/chroma-meas/lib:$LD_LIBRARY_PATH

WORKDIR /meas
# CMD ["bash"]
