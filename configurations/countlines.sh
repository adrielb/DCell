find . -regex '.*\.\(c\|h\|cxx\|cpp\)' -print0 | xargs -0 cat | wc -l

