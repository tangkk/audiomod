# Use CentOS 8 as the base image
FROM centos:8

# Update repository URLs to use vault.centos.org
RUN sed -i 's|mirrorlist.centos.org|vault.centos.org|g' /etc/yum.repos.d/CentOS-*.repo && \
    sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-*.repo

# Install dnf-plugins-core to enable config-manager
RUN dnf install -y dnf-plugins-core && \
    dnf config-manager --set-enabled powertools && \
    dnf groupinstall -y "Development Tools" && \
    dnf install -y epel-release && \
    dnf install -y cmake gcc-c++ glibc-static libstdc++-static boost-devel && \
    dnf clean all

# Set working directory
WORKDIR /app

# Default command
CMD ["bash"]