###### DESCRIPTION
The documentation is generated using Sphinx
     http://www.sphinx-doc.org/en/stable/
together with the matlabdomain extension
     https://pypi.python.org/pypi/sphinxcontrib-matlabdomain
Documentation is generated from Matlab's comments

###### COMPILATION

	sphinx-build -b html source _build
	
	after compiling open the file index.html (in build) and remove the 
	<!-- Local TOC -->
    <div class="local-toc"><ul>
       ...
    </div>
    but let the content that is between these two div (I don't known why there is this problem. It allows to hide subcontents in the TOC)
    
    add style="text-align:justify" in the body "basile":
    <body class="wy-body-for-nav" role="document" style="text-align:justify">
    I did not find a nicer way to do that... 
    
    add <img src="logo.png"  width="300" /> at the following position
    <nav data-toggle="wy-nav-shift" class="wy-nav-side">
    <div class="wy-side-scroll">
       <img src="logo.png"  width="300" />  <<<------- HERE
    <div class="wy-side-nav-search">
    
    set the href of  <a href="#" class="icon icon-home"> GlobalBioIm Library </a> to 
    http://bigwww.epfl.ch/algorithms/globalbioim/

###### HOW TO USE

	- within your comments you can refer to any:
	     - class by      :class:`classname`
	     - attribute by  :attr:`attributename`
	     - method by     :meth:`mathodname`
	
	- You can insert latex within $$ $$ or \\(  \\) (inline) and this will be correctly interpreted.
	  Note that backslashes have to be doubled '\\' (i.e. \int_x --> \\int_x)
	  
	- Let an empty line (i.e. without any %) before the copiright paragraph:
	% See also :class:`LinOp`.

    %     Copyright (C) 2017 ...
    
    
    
	
	
	

