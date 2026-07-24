# common.smk -- shared path helpers used across the rule files.


def get_template_prefix(root, subj_wildcards, template):
    """creates prefix for template files, including subject/session wildcards
    so that DAGs for each subject/session are kept independent.
        e.g.: sub-001/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym"""

    path_entities = bids(root=root, **subj_wildcards).split("/")[
        :-1
    ]  # leave out the file prefix

    path_entities.append(f"tpl-{template}")  # sub-folder name
    path_entities.append(f"tpl-{template}")  # file prefix

    return "/".join(path_entities)
