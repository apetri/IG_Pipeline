#This is a python script to make the plane randomizer file#
from __future__ import print_function,with_statement
from random import seed,randint

#Defaults
number_of_ics = 5
planes_per_snapshot = 9
number_of_snapshots = 60
max_realizations = 15000
mirrot_max = 7
MAX_random_number = 10000000
random_seed = 0

#Parse new values from command line

#Number of ics
try:
	n = int(raw_input("Number of initial conditions to mix (default {0}): ".format(number_of_ics)))
	number_of_ics = n
except ValueError:
	pass

#Planes per snapshot
try:
	n = int(raw_input("How many planes per snapshot (default {0}): ".format(planes_per_snapshot)))
	planes_per_snapshot = n
except ValueError:
	pass

#Number of snapshots
try:
	n = int(raw_input("Number of snapshots (default {0}): ".format(number_of_snapshots)))
	number_of_snapshots = n
except ValueError:
	pass

#Number of ics
try:
	n = int(raw_input("Max number of realizations (including subfields) (default {0}): ".format(max_realizations)))
	max_realizations = n
except ValueError:
	pass

#Mirrot max
try:
	n = int(raw_input("Mirrot max value (default {0}): ".format(mirrot_max)))
	mirrot_max = n
except ValueError:
	pass

#MAX_random_number
try:
	n = int(raw_input("MAX random number (default {0}): ".format(MAX_random_number)))
	MAX_random_number = n
except ValueError:
	pass

#Random seed
try:
	n = int(raw_input("Random seed for generation (default {0}): ".format(random_seed)))
	random_seed = n
except ValueError:
	pass


#Log used values
print("")
print("Mixing {0} initial conditions".format(number_of_ics))
print("5 colums per snapshot x {0} snapshots = your randomizer file will have {1} columns".format(number_of_snapshots,5*number_of_snapshots))
print("This file can handle at most {0} different realizations".format(max_realizations))
print("Mirrot max is {0}".format(mirrot_max))
print("Max random number is {0}".format(MAX_random_number))
print("Initializing with random seed {0}".format(random_seed))
print("")

#Initialize random seed
seed(random_seed)

#Parse filename
filename = raw_input("Where do you want to save the randomizer file? -> ")
print("Writing randomizer info to {0}".format(filename))

#Write info to file
with open(filename,"w") as randomizer:

	for realization in range(max_realizations):
	
		for snapshot in range(number_of_snapshots):
			randomizer.write("{0} {1} {2} {3} {4} ".format(randint(1,number_of_ics),randint(1,planes_per_snapshot),randint(0,mirrot_max),randint(1,MAX_random_number),randint(1,MAX_random_number)))

		randomizer.write("\n")

#Done
print("Done writing {0}".format(filename))
