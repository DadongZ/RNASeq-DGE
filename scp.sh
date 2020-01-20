
key=/mnt/d/github/clouds/aws/keys/me-keypair-ec2-nva.pem

scp -ri ${key} /mnt/d/github/RNASeqDGE/* ubuntu@ec2-3-83-139-154.compute-1.amazonaws.com:/srv/shiny-server/RNASeqDGE/
