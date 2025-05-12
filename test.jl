using PkgTemplates

t = Template(;
           user="hersheytee",
           authors=["Harsh Talathi"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )

t("DetailedBalance")