---
layout: lesson
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---


**Welcome** to this lesson on the fundamental principles of Topological Data Analysis, hereafter called **TDA**. 

TDA is used to analyze large datasets and discover structures that differentiate datasets.

The main objects in TDA are called simplicial complexes and are a generalization of a graph. They are formed by connecting vertices (0-simplices), edges (1-simplices), triangles (2-simplices), and simplices in higher dimensions.

We will begin by exploring the basic concepts of TDA and the computational tools that allow us to use TDA.

In the following chapters, we will delve into the mini Streptococcus database used in the Pangenome Analysis in Prokaryotes lesson to reconstruct a pangenome using simplicial persistence techniques. This will help us identify gene families. In the subsequent lesson, the phenomenon of horizontal gene transfer is studied by searching for 1-holes using two different simplicial complexes.

One of the key highlights of this course is the opportunity to explore how to apply TDA methods to Comparative Genomics analyzes.

Get ready to embark on this exciting journey into the world of TDA.

<!-- this is an html comment -->

{% comment %} This is a comment in Liquid {% endcomment %}

> ## Prerequisites
>
> Before diving into this lesson, it is essential to have a working understanding of the **Bash shell** and the **language Python**. If you are not already familiar with these programming languages, we recommend completing the [Introduction to the Command Line for Pangenomics](https://czirion.github.io/shell-pangenomics/) and [Introduction to Python](https://czirion.github.io/pangenomics-python/).
> 
{: .prereq}

This lesson is the third part of the [Pangenomics Workshop](https://czirion.github.io/pangenomics-workshop/).

{% include links.md %}
