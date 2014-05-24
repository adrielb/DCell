" since env vars not sourced:
let g:petsc_dir="/home/abergman/apps/petsc"
let g:slime_python_ipython = 1

set nocscopeverbose
cscope add cscope.out
"cscope add ~/apps/openmpi-1.8.1/cscope.out
cscope add ~/apps/petsc/cscope.out
set cscopeverbose
set tags+=~/apps/openmpi-1.8.1/tags
"set tags+=~/apps/petsc/CTAGS
set path+=~/apps/petsc
set path+=~/apps/openmpi-1.8.1

map <leader>0 :cd ~/projects/DCell<CR>
map <leader>1 :CtrlP ~/apps/petsc<CR>
map <leader>2 :find local.vimrc<CR>
map <leader>3 :find /home/abergman/tmp/info.log.0<CR>
map <leader>4 :find ./sims/Fibers/Fibers.c<CR>

augroup setDCellDir
  autocmd!
  autocmd InsertLeave * :cd ~/projects/DCell
augroup END

augroup setNomodPETSc
  autocmd!
  autocmd BufEnter ~/apps/petsc/* setlocal nomodifiable 
augroup END

"if exists("g:did_localvimrc")
"  finish
"endif
"let g:did_localvimrc = 1
" stuff to lead only once

"let &efm = "%m line %l in %f,%-G[0]PETSC ERROR:%m" . efm

function! DebugDCellInfo()
 let s:makeprg_save=&makeprg
 let s:efm_save=&efm

 "let &efm ="%+G%.%#ERROR%m,"

 let &efm ="%C%.%#:%m,"
 let &efm.="%+E%.%#ERROR%m,"
 let &efm.="%-Z%.%#MESSAGE,"
 let &efm.="%-G%.%#"

 let &makeprg="cat $PETSC_TMP/info.log.*"
 execute ":silent make"
 execute ":copen"

 let &makeprg=s:makeprg_save
 let &efm=s:efm_save
endfunction

command! D call DebugDCellInfo()
