$machine_cap  = "90"
$machine_cpus = "8"
$machine_name = "ubuntu-EFDC"
$machine_ram  = "4096"

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/trusty64"
 
  config.vm.synced_folder "./", "/home/efdc/", create: true

  config.vm.provider "virtualbox" do |vb|
    vb.customize ["modifyvm", :id, "--name", $machine_name]
    vb.customize ["modifyvm", :id, "--cpus", $machine_cpus]
    vb.customize ["modifyvm", :id, "--cpuexecutioncap", $machine_cap]
    vb.customize ["modifyvm", :id, "--memory", $machine_ram]
  end

  # execute the setup.sh bash script that installs all required dependencies
  config.vm.provision "shell" do |s|
    s.path =  "./setup.sh"
  end
  config.ssh.forward_agent = true
  config.ssh.forward_x11 = true

end
