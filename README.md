# Hello-CFD
Coding practice for me, a beginner of CFD.

# How to use git and push to github
basic command
---
```bash
git status 
git diff
git add .
git commit -m "xxxx"
```
reference:https://liaoxuefeng.com/books/git/remote/add-remote/index.html

# 1. Solve 1D inviscid compressible Euler Equation
This code comes from 胡偶 2011 to understand how to convert the theory into code. 
> * boundary condition: Dummy Cell
> * flux calculation  : AUSM Scheme
> * restruction       : MUSCL 
> * restrictor        :	Van Albada
> * time discrete     :	4-step Runge-Kutta
