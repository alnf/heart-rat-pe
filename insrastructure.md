# Infrastructure configuration

- [Infrastructure configuration](#infrastructure-configuration)
  - [Configuring storage](#configuring-storage)
  - [Configuring services](#configuring-services)
    - [RStudio server](#rstudio-server)


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

## Configuring services

### RStudio server

