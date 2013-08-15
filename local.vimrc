Bundle 'a.vim'
Bundle 'Valloric/YouCompleteMe'

set autowriteall

" 'K' jumps to Petsc Doc
"set keywordprg=ctags
set makeprg=make\ -C\ /home/abergman/Research/DCell\ all

" run make after saving .c file
autocmd BufWritePost *.c make 


