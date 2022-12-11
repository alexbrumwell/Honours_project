# Notes on the cluster
----------------------------------------------------------------------------------------------
## Accessing the cluster
1.	Download and open the university VPN called Cisco AnyConnect
* The server we are using is called “gucsasa.1.cent.gla.ac.uk”
* If I am off-campus then I use the off-campus option
* Then put in my username (i.e. GUID) and password (GUID password)

2.	Use MobaXterm to connect to the remote server (i.e. cluster)
* `> ssh studentprojects@becker.eng.gla.ac.uk`
* Password: 10pointstogryffindor

3.	Access my folder in the cluster:
* `> cd /shared5/Alex`

>Sidenote: To run a function that can run in the background we do:
* `> screen -S <name of task>`
* exit or ctrl + a + d
* To check the running task: `> screen -r <task_name>`
* To close the screen type: `> exit`

## Installing packages
1.	First activate the Conda environment using:
`> export PATH=/home/opt/miniconda2/bin:$PATH`

2.	Then to install the package do:
`> conda create –n trimgalore –c bioconda trim-galore`


## Using software/packages
To use samtools or other software first we need to activate a conda environment using:
```bash
> export PATH=/home/opt/miniconda2/bin:$PATH source activate popgen #activates conda environment
> source activate samtools #activates package
> source deactivate samtools #deactivates package
```
