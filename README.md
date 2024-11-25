# Hello-CFD
I am a beginner of CFD Coding. Here, I store some code for practice and learn. 

# How to use git and push to github
basic command
---
```bash
git status 
git diff
git add .
git commit -m "xxxx"

git quick use https://www.bootcss.com/p/git-guide/
```
reference:https://liaoxuefeng.com/books/git/remote/add-remote/index.html

# 1. Solve 1D inviscid compressible Euler Equation
This code comes from 胡偶 2011 to understand how to convert the theory into code. 
> * boundary condition: Dummy Cell
> * flux calculation  : AUSM Scheme
> * restruction       : MUSCL 
> * restrictor        :	Van Albada
> * time discrete     :	4-step Runge-Kutta

# 2. PINNs: 1D channel flow example
This code comes from dyFluid. [website](http://dyfluid.com/pinn.html)
The environment of OpenFOAM libtorch is needed here.

# 3. CFD-Course-Tsu Homework
The Homework report of university's CFD-Course.
