.. lsst-task-topic:: lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask

######################################
NumberUnassociatedDiaObjectsMetricTask
######################################

``NumberUnassociatedDiaObjectsMetricTask`` computes the number of DIAObjects that are *not* updated when processing data through source association (as the ``ap_association.numUnassociatedDiaObjects`` metric).
It requires task metadata as input.
While the task can operate at image-level or coarser granularity, the current algorithm may double-count objects and should not be run on multiple visits.

.. _lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask-summary:

Processing summary
==================

``NumberUnassociatedDiaObjectsMetricTask`` reads ``AssociationTask`` statistics from task metadata associated with one or more processed images.
It uses these statistics to count the number of unassociated DIAObjects when processing those images.

.. _lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask

.. _lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask-butler:

Butler datasets
===============

Input datasets
--------------

:lsst-config-field:`~lsst.verify.tasks.metadataMetricTask.MetadataMetricConfig.metadata`
    The metadata of the top-level command-line task (e.g., ``ProcessCcdTask``, ``ApPipeTask``) being instrumented.
    Because the metadata produced by each top-level task is a different Butler dataset type, this dataset **must** be explicitly configured when running ``NumberUnassociatedDiaObjectsMetricTask`` or a :lsst-task:`~lsst.verify.gen2tasks.MetricsControllerTask` that contains it.

.. _lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask

.. _lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ap.association.metrics.NumberUnassociatedDiaObjectsMetricTask
