####################################################################################################
####### Steps to Contextualize the phenotype-specific networks and do perturbation simulation ######
####################################################################################################

This is an in-house developed algorithm for making networks phenotype specific referred to as contextualization (Zickenrott et al., Cell Death & Disease, 2016) . In our current study the phenotypes of interest are healthy and Alzheimer's disease phenotype.

This algorithm has many steps involved in making contextualized networks and then identify key genes in the positive and negative circuits of the network and as a last step induce the perturbation and see their effect on the other genes in the network. In the following, I will explain every step one-by-one:


#######################################################
Step 1: Pre-processing the network (Preprocessor.jar):
#######################################################


The preprocessor is able to generate adequate network representations used as an input for the Genetic Algorithm to contextualize the prior knowledge network (PKN) against two phenotypes.

Problem:
Provided a list of genes, Pathway Studio returns only a network of interactions among genes. Some other databases e.g. Metacore, provide 2 files:
	- A file containing information about the nodes of the network
	- a file containing information about the interactions of the network

Now, the nodes represented in MetaCore are not identical to those provided in the initial gene list.
Pre
Example: Input gene CD3g is mapped to 3 different network objects, Cd3g, CD3G and TARP. Thus, interactions might be reported in the form
	- Cd3g —> Gene A
	- CD3G —> Gene B
	- TARP —> Gene C

However, TARP is not represented in our selected genes and we do possibly not have information about this gene. In particular, missing expression values.

NOTE: As in this case we are using Pathway Studio and it doesn't have this problem of node mapping we will generate an other file called nodeMap.txt in which we will paste the *geneList.txt file twice, separated by tab.

Solution:

Usage: java -jar Preprocessor.jar . nodeMap.txt Interactions.txt hmc_genesList.txt

Parameters:
<1>	Path to the directory which contains the input files

Preprocessor.jar now follows four steps to create adequate networks
	1. The reported interactions are extended. All combinations of interactions among mapped gene symbols are built. (File: MappedInteractions.txt)
	2. The resulting mapped interactions might contain genes that are not represented in the list of genes supplied to Pathway Studio or MetaCore. Thus, all of these interactions are filtered out. (File: FilteredMappedInteractions.txt)
	3. Duplicates are removed (File: FilteredMappedInteractionsNoDuplicates.txt)
	4. An adjacency matrix is derived from the resulting interaction list. (adjacency.txt)


########################################################################
Step 2: Differential Network Analysis (DifferentialNetworkAnalysis.jar):
########################################################################

The Differential Network Analysis uses the Genetic Algorithm to filter out the interactions which are not compatible with given gene expression.

Usage: java -jar DifferentialNetworkAnalysis.jar Expression.txt adjacency.txt GAResult.txt 0 true 1000 50 .

Parameters:
<1>	Gene Expression File (tab delimited)
<2>	adjacency matrix
<3>	Output file for genetic algorithm result
<4>	0
<5>	Boolean: true
<6>	population size: 1000
<7>	number of generations: 50

Output:
GAResult.txt
NetworkPhenotype1.txt
NetworkPhenotype2.txt


You have to provide gene expression of 1 phenotype. For example, in this analysis I gave the differential expression profile which represents the AD phenotype (Expression.txt).

The healthy phenotype will be an exact opposite of these expression values.

NetworkPhenotype1.txt contains the interactions which represent the AD phenotype.


##############################################################
Step 3: Common Network Generator (CommonNetworkGenerator.jar):
##############################################################

Common Network Generator generates the set of interactions which are common to both phenotype specific networks.

Usage: java -jar CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt

Parameters:

<1>	Network Phenotype 1
<2>	Network Phenotype 2
<3>	Output file name: CommonNetworkGenerator_Output.txt


##########################################################################
Step 4: Differential Network Generator (DifferentialNetworkGenerator.jar):
##########################################################################

Differential Network Generator generates the set of interactions which are unique to both phenotype specific networks.

Usage: java -jar DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt

Parameters:

<1>	Network Phenotype 1
<2>	Network Phenotype 2
<3>	Output file name: DifferentialNetworkGenerator_Output.txt


##########################################
Step 5: Cmpute Cycles (ComputeCycles.jar):
##########################################

ComputeCycles.jar extracts the circuits of a given network.

Usage: java -jar ComputeCycles.jar CommonNetworkGenerator_Output.txt Expression.txt pos.txt neg.txt

Parameters:
<1>	Interaction List
<2>	Gene Expression File (tab delimited)
<3>	Output file for positive circuits: pos.txt
<4>	Output file for negative circuits: neg.txt


##########################################################
Step 6: Steady State Calculator (SteadyStateCalculator.jar):
##########################################################

SteadyStateCalculator.jar computes the steady state of network starting from a given expression pattern.

Usage:
java -jar SteadyStateCalculator.jar Expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar SteadyStateCalculator.jar Expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt

Parameters:
<1>	Expression file
<2>	Interaction list representing a network
<3>	Either 1 or 2 representing the phenotype the network represents
<4>	Output file

NOTE: Parameter <3> has to be used in the following way. If the input expression is the differential 
expression pattern of one of the two conditions obtained from any criterion and parameter <2> represents 
the phenotype i network, <3> should be i. If the expression pattern is already phenotype specific, 
<3> should be always 1.


##################################################################
Step 6: Perturbagen List Generator (PerturbagenListGenerator.jar):
##################################################################

PerturbagenListGenerator.jar determines genes that might be suitable targets for perturbation 
based on a two-fold decision mechanism. First, nodes residing in circuits of the network are 
collected. As a „measure of importance“ the percentage of circuits each node belongs to is given. 
Nodes with a value closer to 1 should be prioritized since they have the potential to flip more 
genes than the ones with values closer to 0. Second, nodes being differentially regulated are 
included, i.e. that do not share the same set of expressed activators/inhibitors UPON PERTURBATION. 
An integer is provided giving information about the significance of the difference. Negative values 
result from more inhibitors or absence of activators in the second network. Positive values result 
from the same reasoning, replacing inhibitors with activators and vice versa.

Usage: java -jar PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

Parameters:

<1>	List of COMMON positive circuits (as generated by ComputeCycles.jar)
<2>	List of COMMON negative circuits (as generated by ComputeCycles.jar)
<3>	Differential Network (as generated by DifferentialNetworkGenerator.jar)
<4>	Steady State gene expression of the network to be perturbed (as generated by SteadyStateCalculator.jar)
<5>	Output File


#######################################################################
Step 7: Brute Force Perturbations (BruteForcePerturbationsUpdated.jar):
#######################################################################

BruteForcePerturbations.jar computes network response to various kinds of perturbations based on a 
given set of perturbagens to draw from.

Usage: java -jar BruteForcePerturbationsUpdated.jar Expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt

Parameters:
<1>	Expression file
<2>	Interaction list representing the network
<3>	Either 1 or 2 representing the phenotype representation of expression file
<4>	List of Perturbagens (e.g. output of PerturbagenListGenerator.jar)
<5>	Maximum number of genes contained in the perturbation set
<6>	Maximum number of combinations to test for each size of perturbation sets
<7>	Output file

NOTE: Parameter <3> has to be used in the following way. If the input expression is the differential 
expression pattern of one of the two conditions obtained from any criterion and parameter <2> represents 
the phenotype i network, <3> should be i. If the expression pattern is already phenotype specific, 
<3> should be always 1.

If parameter <5> is larger then the number of genes in the perturbagen list, all the genes are taken. 
The maximum number of combinations to take for each size of perturbation sets is bound by the binomial 
coefficient. Thus, if Parameter <6> is 1.000.000 and there are n genes supplied as perturbagens, 
only n sets containing one perturbagen are taken. More formal, if there are n perturbagens then for 
each 1 <= k <= <5>, min(binom(n,k),<6>) combinations will be computed. If <6> is sufficiently large, 
a maximum of pow(2,n) combinations will be tested. (TAKES LONG TIME AND NEEDS A LOT OF MEMORY!! 
NOT RECOMMENDED)
