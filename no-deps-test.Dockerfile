FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8

RUN apt-get update && apt-get install -y \
        ca-certificates \
        curl \
        bash \
        perl-base \
        python3 \
        python-is-python3 \
        git \
        build-essential

WORKDIR /work

CMD ["bash"]
