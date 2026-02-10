import sys
from .__init__ import __version__

def print_logo():
    # Force UTF-8 encoding for stdout
    sys.stdout.reconfigure(encoding='utf-8')
    
    tool_name = "naiveATLAS"
    version = f"{__version__}"
    github_url = "https://github.com/TimothyStephens/naive_atlas"
    owner_email = "ts942@sebs.rutgers.edu"
    
    # Use raw string r'' to handle backslashes
    print("\n\n" + "="*90 + "") # A separator line
    logo = r"""
                    ..
      _ __    __ _  _   _   _   ___     █████╗ ████████╗██╗      █████╗ ███████╗
     | '_ \  / _` || | | | | | / _ \   ██╔══██╗╚══██╔══╝██║     ██╔══██╗██╔════╝
     | | | | |(_| || | | |_| | | __/   ███████║   ██║   ██║     ███████║███████╗
     |_| |_| \__,_||_|  \___/  \___|   ██╔══██║   ██║   ██║     ██╔══██║╚════██║
                                       ██║  ██║   ██║   ███████╗██║  ██║███████║
                                       ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
    """
    print(logo)
    print(f"    {tool_name} (v{version})\n")
    print(f"    GitHub:   {github_url}")
    print(f"    Email:    {owner_email}")
    print("\n\n" + "="*90 + "\n\n") # A separator line

