1. Download the miPIE files from github.com/jrgreen7/miPIE 

2. Download and install miRDeep2 and all dependencies. See https://www.mdc-berlin.de/8551903/en/ for details on installing the miRDeep2 package.

3. Download and install the numpy (http://www.numpy.org/) and scikit-learn (http://scikit-learn.org/stable/) python packages.

4. Within your miRDeep2 installation, replace the miRDeep2.pl and mirdeep_core_algorithm.pl files with those provided in the miPIE package.

5. Run the miRDeep2 pipeline on your NGS data set, as directed by the miRDeep2 instructions.

6. Run the predict.py python script from the base folder of the miPIE installation. The script requires one input, -d, which is the location of the output.mrd file generated in step 5 of these instructions. A typical command line for this script is as follows:

python predict.py -d /home/user/dataset/mirdeep_runs/run_XX_XX_2015_t_XX_XX_XX

Prediction is performed using a pre-generated random forest model which was trained on five NGS data sets.

