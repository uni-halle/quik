# Install Press 2022
sh ./install_cuda_drivers.sh
sh ./init_venv.sh
# Make a RAM-Filesystem for communication with Press 2022
sudo mkdir /media/ramdisk
sudo chmod 777 /media/ramdisk
sudo mount -t tmpfs -o size=20% none /media/ramdisk
# Compile main program
mvn clean package
