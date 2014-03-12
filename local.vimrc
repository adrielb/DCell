" 'K' jumps to Petsc Doc
"set keywordprg=ctags

set sessionoptions=buffers,curdir

" run make after saving .c file
augroup savemake
  autocmd!
  autocmd BufWritePost *.c make
augroup END

cd ~/projects/DCell
set noautochdir
set nocscopeverbose
cscope add cscope.out
cscope add ~/apps/petsc/cscope.out ~/apps/petsc
set cscopeverbose
set tags+=~/apps/petsc/CTAGS

map <leader>1 :find makefile<CR>
map <leader>2 :find local.vimrc<CR>
map <leader>3 :find /home/abergman/tmp/info.log.0<CR>
map <leader>4 :find ./sims/Fibers/Fibers.c<CR>


