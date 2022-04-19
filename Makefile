#
# This is only a very simple makefile for documentation generation
# The real makefile has to be generated through CMake
#

# SOURCE CODE DIRECTORIES LAYOUT
# Set the layout of the source code directory
DOCDIR		= ./Doc/
CODEDOCDIR	= ./CodeDoc/

# DOCUMENTATION CREATION RULES
# Set make rules for documentation generation
#--------------------------------------------------------------------
.PHONY: all
all:
	@echo "#############################################################"
	@echo "#############################################################"
	@echo "#   This Makefile only generates the code documentations.   #"
	@echo "#    - Use 'make doc' to generate the documentation.        #"
	@echo "#    - Use 'make clean' to delete the documentation.        #"
	@echo "#    - Use CMake to generate the real Makefiles.            #"
	@echo "#############################################################"
	@echo "#############################################################"

.PHONY: doc
doc:	
	doxygen
#--------------------------------------------------------------------

#CLEANING RULES
#Set make rules for cleaning source directory
#--------------------------------------------------------------------
.PHONY: clean
clean:
	rm -r -f $(CODEDOCDIR)latex
	rm -r -f $(CODEDOCDIR)html
#--------------------------------------------------------------------
