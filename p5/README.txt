To run cifar10_cnn.py on Bowdoin HPC

cd into the directory
scp cifar10_cnn.py username@bowdoin.edu:~/
ssh to moosehead.bowdoin.edu

hpcsub -l gpu=1 -cmd python cifar10_cnn