.PHONY: example

example:
	rm -f CALYPSO-Bohrium-example.zip
	cp -r example CALYPSO-Bohrium-example
	zip -rq CALYPSO-Bohrium-example.zip CALYPSO-Bohrium-example
	rm -r CALYPSO-Bohrium-example
