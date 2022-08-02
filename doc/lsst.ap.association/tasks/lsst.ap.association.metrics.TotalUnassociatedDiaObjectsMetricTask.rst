.. lsst-task-topic:: lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask

#####################################
TotalUnassociatedDiaObjectsMetricTask
#####################################

``TotalUnassociatedDiaObjectsMetricTask`` computes the number of DIAObjects that have only a single source (as the ``ap_association.totalUnassociatedDiaObjects`` metric).
It requires an alert production database as input, and is meaningful only at dataset-level granularity.

.. _lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask-summary:

Processing summary
==================

``TotalUnassociatedDiaObjectsMetricTask`` queries the database (through `~lsst.dax.apdb.Apdb`) for the number of DIAObjects with exactly one source.

.. _lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask

.. _lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask-butler:

Butler datasets
===============

Input datasets
--------------

:lsst-config-field:`~lsst.verify.tasks.apdbMetricTask.ApdbMetricConfig.dbInfo`
    The Butler dataset from which the database connection can be initialized.
    The type must match the input required by the :lsst-config-field:`~lsst.verify.tasks.apdbMetricTask.ApdbMetricConfig.dbLoader` subtask (default: the top-level science task's config).

.. _lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask

.. _lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ap.association.metrics.TotalUnassociatedDiaObjectsMetricTask
