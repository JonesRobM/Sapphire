# SAPPHIRE
Post - Processing Software
![Sapphire-logos_black](https://user-images.githubusercontent.com/52043020/154812758-c0184aa5-c2d0-4b5c-b23f-5955221d1d79.png)

Sapphire is a multi-faceted platform engineered to facilitate the design, characterisation, and classification of metallic nanoparticles in addition to computing their dynamics and energetics. Integrating these two aspects of nanoparticle simulation is the comprehensive set of analysis techniques offered by Sapphire which utilises pre-exisiting, highly-efficient, pythonic libraries such as ase (atomic simulation environment.

Please, also be aware that Sapphire is still in a perpetual state of development. Should you recover any strange results, contact one of the authors of this document with any information you are able to provide about the problem and a hot-fix will be implemented if necessary. Otherwise, feedback on ease of use is always appreciated.

Installation instructions:

1: Create a python virtual environment

$ python3 -m pyvenv < MyEnv >

2: Source the environment

$ source /path/to/< MyEnv >

3: Link your jupyter notebook to this new environment

$ pip install --user ipykernel

$ python -m ipykernel install --user --name=myenv

Check the file /home/user/.local/share/jupyter/kernels/< MyEnv >/kernel.json

It should look like this

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

{
 "argv": [
  "/path/to/< MyEnv >",
  "-m",
  "ipykernel_launcher",
  "-f",
  "{connection_file}"
 ],
 "display_name": "< MyEnv >",
 "language": "python"
}

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


4: Build Sapphire - its dependencies should be installed for you
Do this from the base directory of Sapphire where you will find the setup.py file

$ python -m pip install --upgrade .

5: Sapphire should now appear in your site packages and will be useable via

$ python

>>import Sapphire

>>from Sapphire.Post_Process import * 

etc...


To get started, visit the Tutorials folder which will guide you through various different utilities of Sapphire 

Please check out our group page to learn more about the exciting work we do at the nanoscale! 

http://balettogroup.org

Authors:

1. Robert M. Jones[^1].

2. Matteo Tibberi

3. Armand Aquier

4. Claudio Zeni

5. Francesca Baletto

[^1]: Robert.M.Jones@kcl.ac.uk
