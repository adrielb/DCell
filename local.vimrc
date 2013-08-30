Bundle 'a.vim'
Bundle 'Valloric/YouCompleteMe'

set autowriteall
autocmd filetype c set cindent

" 'K' jumps to Petsc Doc
"set keywordprg=ctags
set makeprg=make\ -C\ /home/abergman/Research/DCell\ all

" run make after saving .c file
augroup savemake
  autocmd!
  autocmd BufWritePost *.c make
augroup END
