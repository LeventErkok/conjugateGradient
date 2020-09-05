# (c) Copyright Levent Erkok. All rights reserved.
#
# The conjugateGradient library is distributed with the BSD3 license. See the LICENSE file
# in the distribution for details.
SHELL     := /usr/bin/env bash
DEPSRCS   = $(shell find . -name '*.hs' -or -name '*.lhs' -or -name '*.cabal' | grep -v Paths_conjugateGradient.hs)
CABAL     = cabal
TIME      = /usr/bin/time

define mkTags
	@find . -name \*.\*hs | xargs fast-tags
endef

.PHONY: all install sdist clean docs hlint tags

all: install

install: $(DEPSRCS) Makefile
	$(call mkTags)
	@$(CABAL) new-install --lib

test: install
	@echo "*** Starting inline tests.."
	@$(TIME) doctest Math/ConjugateGradient.hs --fast --no-magic -package random

sdist: install
	@(set -o pipefail; $(CABAL) sdist)

veryclean: clean
	@-ghc-pkg unregister conjugateGradient

clean:
	@rm -rf dist

docs:
	cabal new-haddock --haddock-option=--hyperlinked-source --haddock-option=--no-warnings	

release: clean install sdist hlint test docs
	@echo "*** conjugateGradient is ready for release!"

hlint:
	@echo "Running HLint.."
	@hlint Math -i "Use otherwise" -i "Parse error"

ci:
	haskell-ci conjugateGradient.cabal --no-tests --no-benchmarks --no-doctest --no-hlint --email-notifications --no-haddock

tags:
	$(call mkTags)
