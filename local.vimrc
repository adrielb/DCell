"Bundle 'a.vim'
"Bundle 'Valloric/YouCompleteMe'

" 'K' jumps to Petsc Doc
"set keywordprg=ctags
"set makeprg=make\ -C\ /home/abergman/Research/DCell\ all

" run make after saving .c file
augroup savemake
  autocmd!
  autocmd BufWritePost *.c make
augroup END

CDC
set nocscopeverbose
cscope add cscope.out
cscope add ~/apps/petsc/cscope.out ~/apps/petsc
set cscopeverbose
set tags+=~/apps/petsc/CTAGS
