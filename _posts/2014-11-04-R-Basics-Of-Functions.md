---
layout: post
title: "Activity: Working With Functions"
date:   2014-10-29 20:49:52
authors: Adapted from Software Carpentry Materials by Leah A. Wasser
categories: [Using R]
tags : []
description: "This activity will provide the basics of creating a function in R."
code1: 
image:
  feature: textur2_FieldWork.jpg
  credit: Ordway Swisher Biological Station, NEON, thx to Courtney Meier
  creditlink: http://www.neoninc.org
permalink: /Using-R/Working-With-Functions
---
<section id="table-of-contents" class="toc">
  <header>
    <h3 >Contents</h3>
  </header>
<div id="drawer" markdown="1">
*  Auto generated table of contents
{:toc}
</div>
</section><!-- /#table-of-contents -->

# Creating Functions

If we only had one data set to analyze, it might be quick to load the file into a spreadsheet and use that to plot some statistics. 
But what if we have 10 or 12 or more files files to analyze? In addition, that number of files may increase through time as more data are collected. In this lesson, we'll learn how to write a function so that we can repeat several operations with a single command.

## Objectives

* Define a function that takes parameters (input values).
* Return a value from a function.
* ???Test and debug a function.???
* Set default values for function parameters.
* Explain why we should divide programs into small, single-purpose functions.
* Learn how to define a Function

##getting started
Let's start by defining a function `fahr_to_kelvin` that converts temperature values from Fahrenheit to Kelvin:

{% highlight r %}
fahr_to_kelvin <- function(temp) {
    kelvin <- ((temp - 32) * (5/9)) + 273.15
    kelvin
}
{% endhighlight %}

The definition opens with the name of your new function, which is followed by the call to make it a `function` and a parenthesized list of parameter names. The parameters are the input values that the function will use to perform any calculations. In this case, the input might be the temperature value that we wish to convert from fahrenheit to kelvin. You can have as many input parameters as you would like (but too many might be bad style). The body, or implementation, is surrounded by curly braces `{ }`. In many languages, the body of the function - the statements that are executed when it runs - must be indented, typically using 4 spaces. 

HINT: While it is not mandatory to indent your code 4 spaces within a function, it is  strongly recommended as good practice!

When we call the function, the values we pass to it are assigned to those variables so that we can use them inside the function. The last line within the function is what R will evaluate as a returning value. Remember that the last line has to be a command that will print to the screen, and not an object definition, otherwise the function will return nothing - it will work, but will provide no output. For example, let's try running our function. Calling our own function is no different from calling any other function:

{% highlight r %}
fahr_to_kelvin(32)
paste('boiling point of water:', fahr_to_kelvin(212))
{% endhighlight %}

We've successfully called the function that we defined, and we have access to the value that we returned. However, it the function was redefined as follows

{% highlight r %}
fahr_to_kelvin <- function(temp) {
    kelvin <- ((temp - 32) * (5/9)) + 273.15
}
{% endhighlight %}

Now typing

{% highlight r %}
fahr_to_kelvin(32)
{% endhighlight %}

Will return nothing.

--> In Python lessons, now would come debugging…
==============================


## Composing Functions

Now that we've seen how to turn Fahrenheit into Kelvin, it's easy to turn Kelvin into Celsius:

{% highlight r %}
kelvin_to_celsius <- function(temp) {
    Celsius <- temp - 273.15
	Celsius
}

paste('absolute zero in Celsius:', kelvin_to_celsius(0))
{% endhighlight %}

What about converting Fahrenheit to Celsius? We could write out the formula, but we don't need to. Instead, we can compose the two functions we have already created:

{% highlight r %}
fahr_to_celsius <- function(temp) {
	temp_k <- fahr_to_kelvin(temp)
	result <- kelvin_to_celsius(temp_k)
    result
}

paste('freezing point of water in Celsius:', fahr_to_celsius(32.0))
{% endhighlight %}

This is our first taste of how larger programs are built: we define basic operations, then combine them in ever-large chunks to get the effect we want. 
Real-life functions will usually be larger than the ones shown here—typically half a dozen to a few dozen lines—but they shouldn't ever be much longer than that, or the next person who reads it won't be able to understand what's going on. __Modular programming__