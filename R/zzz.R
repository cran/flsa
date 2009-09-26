########################################################
###
### Define .First.lib and .Last.lib functions for the package
###
########################################################

.First.lib = function(libname, pkgname)
{
    library.dynam("flsa", package=pkgname)
}

.Last.lib = function(libpath)
{
    library.dynam.unload("flsa", libpath)
}

