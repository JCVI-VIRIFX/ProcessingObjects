#
#  Makefile for PCR resequencing data analysis applications 
#           Processing Objects perl modules
#  $Id$
#

CURDIR      = .
INSTALLDIR  = /usr/local/reseq/data_analysis/ProcessingObjects
INSBINDIR   = $(INSTALLDIR)/bin
INSLIBDIR   = $(INSTALLDIR)/lib
INSETCDIR   = $(INSTALLDIR)/etc
INSDOCDIR   = $(INSTALLDIR)/doc

install:
	@mkdir -p $(INSBINDIR)
	@cp Metrics.pm $(INSBINDIR)/
