all:
	(cd serial/laspack; $(MAKE));
	(cd serial; $(MAKE));

rm:
	(cd serial/laspack; $(MAKE) rm);
	(cd serial; $(MAKE) rm);
