R=R
INPUT= $(shell find . -maxdepth 1 -type f -name "*.Rmd")
MD_FILES= $(patsubst %.Rmd,%.md,$(INPUT))
HTML_FILES= $(patsubst %.Rmd,%.html,$(INPUT))

all: setup $(MD_FILES)

html: all $(HTML_FILES)

setup:
	mkdir -p gmt
	mkdir -p compare_signature_files

clean:
	rm -fv *.html
	find . -maxdepth 1 -type f -name "*.md" | grep -v README | xargs rm -fv
	rm -rfv *_files
	rm -rfv fig/*
	rm -rfv gmt/*

wipe: clean
	rm -rfv *_cache

%.md: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"

%.html: %.md
	pandoc -s --mathjax --number-sections -o $@ $<
