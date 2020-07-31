# Makefile for qutip-bo

redb = "\033[1;31m"$1"\033[0m"
greenb = "\033[1;32m"$1"\033[0m"
green = "\033[0;32m"$1"\033[0m"
yellow = "\033[0;33m"$1"\033[0m"
whiteb = "\033[1;37m"$1"\033[0m"
DOCKERHUB_ORGANISATION=julianiacoponi
VERSION?=$(shell ./inc-version)

help:
	@echo $(VERSION)
	@echo $(call greenb,"Please use make <target>' where <target> is one of:")
	@echo $(call redb," NOTHING IMPLEMENTED PLEASE DON'T USE!")" (please...)"

pull-%:
	docker pull $(DOCKERHUB_ORGANISATION)/$*
export_version:
	@echo "__version__ = '$(VERSION)'" > thetis/__init__.py
	@echo $(VERSION) > .VERSION

qutip-bo: export_version
	docker build --tag=${DOCKERHUB_ORGANISATION}/$@:$(VERSION) -f dockerfiles/Dockerfile.qutip-bo .
	docker tag ${DOCKERHUB_ORGANISATION}/$@:${VERSION} $@
	docker tag ${DOCKERHUB_ORGANISATION}/$@:${VERSION} ${DOCKERHUB_ORGANISATION}/$@

push-%:
	docker tag $(DOCKERHUB_ORGANISATION)/$*:$(VERSION) $(DOCKERHUB_ORGANISATION)/$*:latest
	docker push $(DOCKERHUB_ORGANISATION)/$*:$(VERSION)
	docker push $(DOCKERHUB_ORGANISATION)/$*:latest

%.docker_image: %
	docker tag $(DOCKERHUB_ORGANISATION)/$*:$(VERSION) $*:latest
	docker save --output=$@ $*:latest

dev-image: export_version
	@echo  $(call green," Launch dev image")
	@if [ ! -f key.pem ] ; then \
		openssl req -x509 -newkey rsa:2048 -keyout key.pem -out cert.pem -nodes -subj "/C=GB/L=London/O=thetis/CN=thetis.example.com"; \
		fi
	@if [ $(latest) ] ; then make pull-thetis; fi
	docker build --tag=dev -f configurations/bedfont-testenv/Dockerfile.dev-image .
	@# touch the file so that docker doesn't make a directory for it
	@touch $(shell pwd)/.dev-image-history
ifdef BISTON_PATH
	docker run --privileged -it --net="host" --rm=true -v $(shell pwd):/app/ -v $(BISTON_PATH):/biston/ -v $(shell pwd)/core:/core -v /mnt:/mnt -v /lib/modules:/lib/modules -v $(shell pwd)/.dev-image-history:/root/.bash_history --entrypoint=/bin/bash thetis_dev
else
	docker run --privileged -it --net="host" --rm=true -v $(shell pwd):/app/ -v $(shell pwd)/core:/core -v /mnt:/mnt -v /lib/modules:/lib/modules -v $(shell pwd)/.dev-image-history:/root/.bash_history --entrypoint=/bin/bash thetis_dev
endif
	@rm -f std* ptmx kcore null zero tty urandom full console
	@rm -rf fd
	@rm -ff fuse
