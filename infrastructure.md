# Infrastructure configuration

- [Infrastructure configuration](#infrastructure-configuration)
  - [Configuring storage](#configuring-storage)
  - [SSH connection](#ssh-connection)
  - [Configuring services](#configuring-services)
    - [RStudio server](#rstudio-server)
    - [JupyterHub](#jupyterhub)


The work is done using [denbi-cloud.bioquant.uni-heidelberg.de](https://denbi-cloud.bioquant.uni-heidelberg.de) resources.

## Configuring storage

Volume size is limited by 2 TB, so we will use it for scrath area and the rest will be stored on NFS share. First, we need to create NFS and volumes through OpenStack dashboard. When creating share it is important to add IP access rule and specify floating IP of the running VM.

Second, we need to configure [volumes](https://cloud.denbi.de/wiki/simple_vm/volumes/) and [NFS share](https://cloud.denbi.de/wiki/Compute_Center/Heidelberg/).

Check devices:

```bash
lsblk -o NAME,SIZE,MOUNTPOINT,FSTYPE,TYPE | egrep -v "^loop"
```

Format the volume:

```bash
mkfs.ext4 /dev/vdb
```

Mount volume:

```bash
sudo mkdir -p /mnt/scratch
sudo mount /dev/vdb /mnt/scratch
sudo chown -R ubuntu:ubuntu /mnt/scratch
```

Mount NFS share (check complete mount path by cliking on share in the dashboard):

```bash
sudo apt install nfs-common
sudo mkdir -p /mnt/data
sudo mount -o vers=4.0 isiloncl1-487.denbi.bioquant.uni-heidelberg.de:/ifs/denbi/manila-prod/YOUR-SHARE /mnt/data
sudo chown -R ubuntu:ubuntu /mnt/data
```

Add following lines to fstab:

```bash
UUID=uuid_of_your_volume    /mnt/scratch    auto    defaults    0   2
isiloncl1-487.denbi.bioquant.uni-heidelberg.de:/ifs/denbi/manila-prod/YOUR-SHARE    /mnt/data   nfs vers=4.0    0   0
```

Now even after the VM reboot all storages will be automatically mounted.

## SSH connection

In order to connect to the VM directly, without the need to go through the jump host first, one needs to [configure ProxyJump](https://cloud.denbi.de/wiki/Compute_Center/Heidelberg/#connecting-to-your-vms-directly) in the ssh config.

## Configuring services

### RStudio server

The overall deNBI manual is [here](https://cloud.denbi.de/wiki/Tutorials/RStudio_Server/). But we will also follow official RStudio Server installation [manual](https://posit.co/download/rstudio-server/).

Do not forget to set up the password for the default Ubuntu user. After we configure ProxyJump it is easy to ssh port-forward to use RStudio Server.

```bash
ssh ubuntu@ip-address -L 8787:localhost:8787
```

### JupyterHub

JupyterHub [installation](https://jupyterhub.readthedocs.io/en/stable/quickstart.html) is quite straightforward. We will install it inside conda environment to avoid versions conflicts. I've also had to add additional command:

```bash
pip install jupyter
```

Port forwarding is as usual:

```bash
ssh ubuntu@ip-address -L 8000:localhost:8000
```
